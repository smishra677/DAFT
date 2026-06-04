import copy
from collections import Counter
import sys
import pandas as pd
import time
import pprint
import re
import argparse
import os
import pickle
import numpy as np
from concurrent.futures import ThreadPoolExecutor, TimeoutError

from rich.console import Console, Group
from rich.live import Live
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
from rich.table import Table

import warnings

warnings.filterwarnings(
    "ignore",
    message="Downcasting object dtype arrays on .*fillna.*",
    category=FutureWarning,
)

warnings.filterwarnings(
    "ignore",
    message="DataFrameGroupBy.apply operated on the grouping columns.*",
    category=DeprecationWarning,
)


def parse1():
    parser = argparse.ArgumentParser(description="DAFT Tranform")

    parser.add_argument('--sp', type=str, help="Species tree")
    parser.add_argument('--path', type=str, default="./DAFT_utils", help="Path to DAFT_utils")
    #parser.add_argument('--verbose', type=int, default=1, help="Verbose mode. 1 = yes, 0 = no")
    parser.add_argument(
        '--lineages',
        type=lambda s: s.split('/'),
        help="Test lineage separated  by '/' (e.g. l1/l2)"
    )

    parser.add_argument('--gt_list', nargs='+', help="List of gene trees (each a Newick string)")
    parser.add_argument('--gt_list_index', nargs='+', help="List of gene trees index for caching")

    parser.add_argument(
        '--verbose',
        type=int,
        default=0,
        help="Show DAFT/djiNNI progress dashboard. 1 = yes, 0 = no"
    )

    parser.add_argument('--output', type=str, help="Name of output file")
    parser.add_argument('--cache_hash', type=str, default="no_hash", help="Input hash used for djiNNI cache")
    parser.add_argument('--random_seed', type=int, default=42, help="Random seed for djiNNI")
    args = parser.parse_args()
    return args


parser = parse1()

path = parser.path

sys.path.append(path)
sys.path.append(path + "/reconcILS")

from reconcILS import *
from utils_reconcILS import *
from DAFT_essential import *


sp_string = parser.sp
lineages = parser.lineages
output = parser.output

lineage1 = lineages[0]
lineage2 = lineages[1]


lis = []

gt_list = parser.gt_list
gt_list_index = parser.gt_list_index
#show_progress = parser.progress == 1
verbose = parser.verbose == 1


djiNNI_hash=parser.cache_hash
random_seed =parser.random_seed
np.random.seed(random_seed)


essential = daft_essential()
reco = reconcils()
red = readWrite.readWrite()



def print_extension(message,out=output.split('_')[:-2]):
   
    log_path = '_'.join(out) + "_djiNNI_log.txt"

    with open(log_path, "a", encoding="utf-8") as log:
        log.write(message+'\n')
        
   
       
def format_seconds(seconds):
    seconds = int(max(seconds, 0))
    h = seconds // 3600
    m = (seconds % 3600) // 60
    s = seconds % 60
    return f"{h:02d}:{m:02d}:{s:02d}"


def get_progress_checkpoints(total):
    fixed = [10, 20, 50, 100, 250, 500, 750, 1000]
    percentages = [
        int(total * p)
        for p in [0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00]
    ]
    checkpoints = set(fixed + percentages + [total])
    checkpoints = {x for x in checkpoints if 0 < x <= total}
    return checkpoints


class DAFTProgress:
    def __init__(self, total_trees, lineage1, lineage2,verbose):
        self.total_trees = total_trees
        self.lineage1 = lineage1
        self.lineage2 = lineage2
        self.start_wall_time = time.perf_counter()
        self.verbose=verbose

        self.console = None
        self.progress = None
        self.task_id = None
        self.live = None

        # Runtime estimate excludes cache hits.
        self.uncached_time = 0.0
        self.uncached_count = 0
        self.cached_count = 0
        self.completed = 0
        self.closed = False
        self.checkpoints = get_progress_checkpoints(total_trees)
        
        if verbose:
                self.console = Console()
                self.progress = Progress(
                SpinnerColumn(),
                TextColumn("[bold cyan]Reconciling gene trees"),
                BarColumn(),
                TaskProgressColumn(),
                TextColumn("[bold]{task.completed}/{task.total}"),
                TimeElapsedColumn(),
                TimeRemainingColumn(),
                expand=True,
                )
                self.task_id = self.progress.add_task("DAFT/djiNNI", total=total_trees)
                self.live = Live(
                self._render(),
                console=self.console,
                refresh_per_second=4,
                transient=False,
                )
                self.live.start()

    #def record_cached(self):
    #self.cached_count += 1

    def record_uncached(self, seconds):
        self.uncached_count += 1
        self.uncached_time += seconds

    def _estimate_values(self):
        elapsed_wall = time.perf_counter() - self.start_wall_time

        completed = int(self.progress.tasks[self.task_id].completed)
        remaining_total = self.total_trees - completed

        # Forecast 1: actual run speed, including cache
        if completed > 0:
                avg_wall = elapsed_wall / completed
                wall_remaining = avg_wall * remaining_total
                wall_total = avg_wall * self.total_trees

                wall_avg_text = f"{avg_wall:.4f} sec/tree"
                wall_remaining_text = format_seconds(wall_remaining)
                wall_total_text = format_seconds(wall_total)
        else:
                wall_avg_text = "calculating..."
                wall_remaining_text = "calculating..."
                wall_total_text = "calculating..."

        # Forecast 2: uncached reconciliation speed only
        if self.uncached_count > 0:
                avg_uncached = self.uncached_time / self.uncached_count
                uncached_remaining = avg_uncached * remaining_total
                uncached_total = avg_uncached * self.total_trees

                uncached_avg_text = f"{avg_uncached:.4f} sec/tree"
                uncached_remaining_text = format_seconds(uncached_remaining)
                uncached_total_text = format_seconds(uncached_total)
        else:
                uncached_avg_text = "no uncached trees yet"
                uncached_remaining_text = "no uncached estimate"
                uncached_total_text = "no uncached estimate"

        return {
                "elapsed_wall": elapsed_wall,

                "wall_avg_text": wall_avg_text,
                "wall_remaining_text": wall_remaining_text,
                "wall_total_text": wall_total_text,

                "uncached_avg_text": uncached_avg_text,
                "uncached_remaining_text": uncached_remaining_text,
                "uncached_total_text": uncached_total_text,
        }

    def _render(self):
                values = self._estimate_values()

                table = Table.grid(padding=(0, 2))
                table.add_column(justify="left", style="bold")
                table.add_column(justify="right")

                completed = int(self.progress.tasks[self.task_id].completed)

                table.add_row("Significant pair", f"{self.lineage1} vs {self.lineage2}")
                table.add_row("Total gene trees detected", str(self.total_trees))
                table.add_row("Completed gene trees", f"{completed}/{self.total_trees}")
                table.add_row("Elapsed wall time", format_seconds(values["elapsed_wall"]))

                table.add_row("", "")
                table.add_row("avg time/tree", values["wall_avg_text"])
                table.add_row("Remaining time", values["wall_remaining_text"])
                table.add_row("Total runtime", values["wall_total_text"])

                
                return Panel.fit(
                        Group(self.progress, table),
                        title="[bold green]DAFT Direction / djiNNI[/bold green]",
                        border_style="red",
                        padding=(1, 2),
                )

    def update(self, completed):
        self.completed = completed

        if not self.verbose:
            return

        self.progress.update(self.task_id, completed=completed)
        self.live.update(self._render())
        if completed >= self.total_trees:
            self.close()

    def close(self):
        if not self.verbose:
            return

        if self.closed:
            return
        self.closed = True

        if self.live is not None:
            self.live.stop()

        elapsed_wall = time.perf_counter() - self.start_wall_time
        print()
        print("=" * 80)
        print("DAFT Direction / djiNNI completed")
        print("=" * 80)
        print(f"Total gene trees:          {self.total_trees}")
        print(f"Reconciled from scratch:   {self.uncached_count}")
        print(f"Loaded from cache:         {self.cached_count}")
        print(f"Total elapsed wall time:   {format_seconds(elapsed_wall)}")
        print("=" * 80)
          

def nni_search(cur_address,sibling):
        stack =[cur_address]
        if type(sibling)!=set:
               sibling={sibling}
        
        count=0
        while stack:
                current_node =stack.pop()
                
                if current_node:
                        #print('===============',sibling, current_node.taxa,current_node.parent)
                        if current_node.parent:
                                current_taxa=current_node.parent.taxa
                                if type(current_node.parent.taxa)!=set:
                                        current_taxa={current_taxa}
                                if len(sibling.intersection(current_taxa))>0:
                                        another_side= [chil for chil in (current_node.parent.children) if chil!=current_node][0]
                                        if another_side.isLeaf:
                                                return count


                                else:   
                                        #print(2222222222222222)
                                        stack.append(current_node.parent)
                                        count=count+1
                        else:   
                                #print(55555555555555555)
                                another_side= [chil for chil in (current_node.children) if len(sibling.intersection(chil.taxa))>0][0]
                                stack.append(another_side)
        
        stack_2=[another_side]
        while stack_2:
                current_node =stack_2.pop()
                #rint()
                if current_node:
                        #print('xxxxxxxxxxxxxxxxxx',sibling, current_node.taxa,)
                        if current_node.leftChild:
                                current_taxa_left=current_node.leftChild.taxa
                                if type(current_node.leftChild.taxa)!=set:
                                                current_taxa_left={current_taxa_left}
                                if len(sibling.intersection(current_taxa_left))>0:
                                        if current_node.leftChild.isLeaf:
                                                return count
                                        else:
                                                stack_2.append(current_node.leftChild)
                                                count=count+1
                        
                        if current_node.rightChild:
                                current_taxa_right=current_node.rightChild.taxa
                                if type(current_node.rightChild.taxa)!=set:
                                                current_taxa_right={current_taxa_right}
                                if len(sibling.intersection(current_taxa_right))>0:
                                        if current_node.rightChild.isLeaf:
                                                return count
                                        else:
                                                stack_2.append(current_node.rightChild)
                                                count=count+1                        
                        
                
        return count

def nni_possible(taxa,species_tree,gene_tree):
        
        reco =reconcils()
        #red= readWrite.readWrite()

        sibling=reco.get_current_sister(species_tree,taxa)

        #print(species_tree.to_newick(),sibling)
        cur_address =essential.current_address(taxa,gene_tree)
        #print(sibling , cur_address,sibling and cur_address)

        if sibling and cur_address:
                        return nni_search(cur_address,sibling),cur_address
                        #return 1,cur_address
        else:
                        return 0,cur_address


def initialize_counter(lineage1, lineage2):
    lineage_stats = {
        f"{lineage1} > {lineage2}": 0,
        f"{lineage1} < {lineage2}": 0,
        f"{lineage1} == {lineage2}": 0,
    }

    return lineage_stats


def querry_lineage(sorted_dict,lineage1,lineage2,lineage_stats):
        count1=0
        count2=0
        for  lineage in sorted_dict:
                if essential.isequal(lineage1,lineage):
                    count1=sorted_dict[lineage]
                if essential.isequal(lineage2,lineage):
                    count2=sorted_dict[lineage]
        if count1>count2:
                #print(123)
                lineage_stats[f"{lineage1} > {lineage2}"] += 1  
                #print(lineage_stats)
                return lineage1,lineage_stats
        elif count1<count2:
                #print(2323)
                lineage_stats[f"{lineage1} < {lineage2}"] += 1
                #print(lineage_stats)
                return lineage2,lineage_stats
        else:   
                #print('rando')
                #np.random.seed(42)  
                lineage_stats[f"{lineage1} == {lineage2}"] += 1
                return np.random.choice([lineage2,lineage1]),lineage_stats
            

def  get_nni(taxa,test_dic):
        #print(test_dic)
        for tr in test_dic:
                if essential.isequal(taxa.to_newick(),tr):
                        if test_dic[tr]>0:
                                return test_dic[tr]
                        else:
                                return 0

def find_minimal_set(key_dict): 
    
    
    return [list(key_dict.keys())]
       
       

def max_round(value_tracker,moves):
        red= readWrite.readWrite()
        pool={}
        max_=0
        #value_tracker= tracker.values()
        #print(value_tracker)
        if type(moves)!=list:
               moves=[moves]
        for tree in moves:
                
                tree= red.parse(tree)
                for li in value_tracker:
                        for k in li[:-2]:
                                if essential.isequal(tree.to_newick(),k):
                                        if int(li[-1])>max_:
                                                max_=int(li[-1])

                pool[tree.to_newick()]=max_

        #print(pool)
        return pool,max(pool.values())

def find_all_lineage(tree,sp):
        ret=[]
        stack=[tree]

        while stack:
                current_node= stack.pop()
                if current_node:
                        if essential.current_address(current_node.taxa,sp):
                                ret.append(current_node.to_newick())
                                
                          
                        if current_node.leftChild:
                                stack.append(current_node.leftChild)
                        if current_node.rightChild:
                                stack.append(current_node.rightChild)
        return ret
        
                                
              

       

def find_moved(tree,sp,gene_tree,test_dic):
        
        reco =reconcils()
        red= readWrite.readWrite()
        nni_orignal= get_nni(tree,test_dic)
        
        taxa= tree.taxa
        return_list=[]
        for taxa_ in find_all_lineage(tree,sp):
                #print(taxa_)
                tree_ =red.parse(taxa_)
                taxa_=tree_.taxa
                cur_address =essential.current_address(taxa_,gene_tree)   
                #print(cur_address)
                if cur_address:
                        nni=get_nni(tree_,test_dic)
                        if nni is None:
                                nni=0
                        nni+= nni_orignal
                        return_list.append([cur_address.to_newick(),1,nni,cur_address])
                else:
                       return_list.append([None,0,0,cur_address])

        

        return return_list

def not_in(newick,visited_list):
        for topo in visited_list:
              if essential.isequal(newick,topo):
                     return False
        return True
              

def merge_moving(dic,obj):
        dic_obj ={obj.parse(key):value for key,value in dic.items()}
        
        for key,value in dic_obj.items():
               key.label_internal()
        
        visited=[]
        for key1,value1 in dic_obj.items():
                for key2,value2 in dic_obj.items():
                        if key1!=key2:
                                if essential.current_address(key2.taxa,key1) and key2 not in visited and key1 not in visited:

                                       dic_obj[key1]+=dic_obj[key2]
                                       visited.append(key2)
                visited.append(key1)

        
        return {key.to_newick():value for key,value in dic_obj.items()}
                                       
def account_trios(tracker):
        visited=[]
        value_list1=[]
        value_list =list(tracker.values())
        
        for val in value_list:
                if len(val)>1:
                    for val2 in val:
                        value_list1.append(val2)
                else:
                       value_list1.append(val[0])

        value_list=value_list1
        v=[]
        for val in value_list:
                v+=val[:-2]
        value_list=v

        dic= dict(Counter(value_list))
        #new_dic={}
        visited=[]
        to_delete=[]
        dic_keys=list(dic.keys())
        for index1 in range(len(dic_keys)):
                for index2 in range(index1+1,len(dic_keys)):
                        key1=dic_keys[index1]
                        key2=dic_keys[index2]
                        com_key=sorted([index1,index2]) 
                        #print(key1,key2)
                        if key1!=key2 and essential.isequal(key1,key2) and com_key not in visited :
                                visited.append(com_key)
                                dic[key1]= dic[key1]+dic[key2]
                                to_delete.append(key2)


        
        for delete in to_delete:
                del dic[delete]
                
        return visited,value_list1,value_list,dic
                                       

def tabulate(tracker,sp,gt,test_dic,lineage1,lineage2,lineage_stats):
        red= readWrite.readWrite()
        

        #visited_list=[]
        visited,value_list1,value_list,dic=account_trios(tracker)
        

        
    
        dic=dict(sorted(dic.items(),key=lambda item: item[1],reverse=True))
        dic =merge_moving(dic,red)
        #print(dic)
        
        moving_taxas,lineage_stats=querry_lineage(dic,lineage1,lineage2,lineage_stats)
        #print('==>',moving_taxas,lineage1,lineage2)
        overall_return={}
        moving_taxas=[moving_taxas]
        path=0
        #for move_ in moving_taxas:
        #tree_= red.parse(move_)

        #tree_.label_internal()
        #val_nni =get_nni(tree_,test_dic)
        #print('xxxx>',tree_.taxa,sp,gt.to_newick())
        #flag_nni, cur_address =nni_possible(tree_.taxa,sp,gt)                                        
        #loop_= [[move_,flag_nni,val_nni,cur_address]]

        #for move_,flag_nni,val_nni,cur_address in loop_:
        #if cur_address:
        #if  not_in(cur_address.parent.to_newick(),visited_list):
        #print(moving_taxas,lineage2)
        if moving_taxas[0]==lineage2:
                #sibling =[chilee for chilee in cur_address.parent.children if chilee!=cur_address ][0]
                overall_return[path]={'Sibling':lineage1,'What_moved':lineage2,'NNI':1}
                #visited_list.append(cur_address.parent.to_newick())
                
        elif moving_taxas[0]==lineage1:
                overall_return[path]={'Sibling':lineage2,'What_moved':lineage1,'NNI':1}
                                                
                             
                             
              

        return overall_return,lineage_stats
        

def render_lineage_stats(lineage_stats,verbose):
        if verbose:
                print("\nLineage comparison counts:")
        print_extension("\nLineage comparison counts:")
        for key, value in lineage_stats.items():
                if verbose:
                        print(f"{key}: {value}")
                print_extension(f"{key}: {value}")
        if verbose:
                print('##'*40)
        print_extension('##'*40)
        


      

def _parse_list(s: str):
    return [x for x in re.split(r'[,\s/]+', s.strip()) if x]


write_intro={'idx':[],'Replicate':[],'Path':[],'From_Where_moved':[],'Sibling':[],'What_moved':[],'NNI':[]}
write_intro1={'Moves':[],'Replicate':[],'NNI':[]}


# Cache reconciliation output by gene tree index.
# This assumes the same gene tree file and same gene tree order are used across runs.
CACHE_DIR = "djiNNI_cache"
CACHE_RUN_DIR = os.path.join(CACHE_DIR, djiNNI_hash)
os.makedirs(CACHE_RUN_DIR, exist_ok=True)


print_extension("="*40)

print_extension(f"Comparision lineages {lineage1} and {lineage2}")



if verbose:
        print("djiNNI cache hash:", djiNNI_hash)
        print("djiNNI cache directory:", CACHE_RUN_DIR)
        print_extension(f"djiNNI cache hash: {djiNNI_hash}")
        print_extension(f"djiNNI cache directory: {CACHE_RUN_DIR}")

print_extension("="*40)   

 
def get_reconciliation_cache_path(gene_tree_index):
        return os.path.join(CACHE_RUN_DIR, f"gene_tree_{gene_tree_index}.pkl")

lineage_stats= initialize_counter(lineage1,lineage2)

#xit()
#exit()
#exit()
sp_large= red.parse_bio(sp_string)
sp_large.label_internal()
#red.write_introgression(sp_large)
done=[]



progress = None

progress = DAFTProgress(len(gt_list), lineage1, lineage2,verbose)

for tree_number,gene_tree in enumerate(gt_list):
        completed = tree_number + 1
        k=gt_list_index[tree_number]
        #print(k,gene_tree,lineage1, lineage2)
        #if not show_progress:
        #print(k,gene_tree,' '*3, output)
        
        #print('gene tree :',gene_tree , 'Test_Lineage :',lineage1 , lineage2)

        gene_tree_string = gene_tree.replace('e-', '0')
        #print(gene_tree_string)

        tr= red.parse_bio(gene_tree_string)
        #print(tr.to_newick())
        
                
        sp= red.parse_bio(sp_string)
        lis.append(tr.to_newick())


        tr=tr.parse(gene_tree_string)
        sp=sp.parse(sp_string)

        sp_copy= copy.deepcopy(sp)

        start_time1= time.time()

        reco.gene_tree=copy.deepcopy(tr)
        tr.order_gene(sp)
        tr.label_internal()
        sp.label_internal()
        tr.map_gene(sp)
        reco.setCost(sp)
        sp.isRoot=True
        tr.isRoot=True
        reco.introgression=True
        reco.L_cost=  20000000000000
        reco.D_cost=  20000000000000
        sp_copy.isRoot=True

        
        def call_reconcILS_function():
                return reco.iterative_reconcILS(tr,sp,sp_copy,sp,[],{},{})
        
        reco.gene_tree= copy.deepcopy(tr)
        species_edge_list=reco.get_edges(sp)

        if len(gene_tree_string)<0:
                reco.reconcILS(tr,sp,sp_copy,sp,{})
                
        
                li =red.sp_event(sp,[])
        else:
                #print('Using Iterative Function')
                num_threads = 40
                timeout_duration = 1800
                reconciliation_cache_path = get_reconciliation_cache_path(k)

                if os.path.exists(reconciliation_cache_path):
                        #if not show_progress:
                        if verbose:
                                print(f"Loading cached reconciliation for gene tree index {k}")
                        print_extension(f"Loading cached reconciliation for gene tree index {k}")
                        with open(reconciliation_cache_path, "rb") as cache_file:
                                cached_reconciliation = pickle.load(cache_file)

                        li = cached_reconciliation["li"]
                        introgression = cached_reconciliation["introgression"]
                        tracker = cached_reconciliation["tracker"]
                        test_dic = cached_reconciliation["test_dic"]

                        #if progress is not None:
                        #progress.record_cached()
                else:
                        uncached_start = time.perf_counter()

                        with ThreadPoolExecutor(max_workers=num_threads) as executor:
                                future = executor.submit(call_reconcILS_function)

                        try:
                                li , introgression, tracker,test_dic= future.result(timeout=timeout_duration)
                        except TimeoutError:
                                print(f"The function call timed out after {timeout_duration} seconds.")
                                exit()
                                li=[]

                        uncached_elapsed = time.perf_counter() - uncached_start

                        if progress is not None:
                                progress.record_uncached(uncached_elapsed)

                        with open(reconciliation_cache_path, "wb") as cache_file:
                                pickle.dump({
                                        "li": li,
                                        "introgression": introgression,
                                        "tracker": tracker,
                                        "test_dic": test_dic
                                }, cache_file)

                #if progress is not None:
                progress.update(completed)
                
        

                #print('-----------------------------------------------------Done With reconcils-----------------------------------------------------------------------------')
                #import pprint
                #pprint.pprint(li)
                #print('123'*100)
                #pprint.pprint(tracker)
                
                if len(tracker)>=1:
                        table,lineage_stats=tabulate(tracker,sp,tr,test_dic,lineage1,lineage2,lineage_stats)
                        added=0
                        #exit()
                        #print(table)
                        for path in table:
                                write_intro['idx']+=[k]
                                write_intro['Replicate']+=[reco.gene_tree.to_newick()]
                                write_intro['Path']+=[path]
                                write_intro['NNI']+=[1]#table[path]['NNI']]
                                write_intro['From_Where_moved']+=["("+lineage2[:-1]+","+lineage1[:-1]+")"]
                                write_intro['What_moved']+=[table[path]['What_moved']]
                                write_intro['Sibling']+=[table[path]['Sibling']]
                else:   
                             
                        #print(gene_tree_string,tracker)
                        #print('We are here')
                        #exit() 
                        '''  
                        write_intro['Replicate']+=[reco.gene_tree.to_newick()]
                        write_intro['Path']+=[0]
                        write_intro['NNI']+=[0]
                        write_intro['From_Where_moved']+=['N/A']
                        write_intro['What_moved']+=['N/A']
                        write_intro['idx']+=[k]
                        write_intro['Sibling']+=['N/A']
                        '''
                        continue  
                                
                intro_df= pd.DataFrame(write_intro)
                intro_df.to_csv('djiNNI_'+output+'.csv', index=False)
                
#if verbose:
render_lineage_stats(lineage_stats,verbose)