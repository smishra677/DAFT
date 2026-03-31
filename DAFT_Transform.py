import copy
from collections import Counter
import sys
sys.path.append("./DAFT_utils")
sys.path.append("./DAFT_utils/reconcILS")
from reconcILS import *
from utils_reconcILS import *
import pandas as pd 
import time
from DAFT_essential import *
import pprint
import re
import argparse
import numpy as np
from concurrent.futures import ThreadPoolExecutor, TimeoutError


          
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



def querry_lineage(sorted_dict,lineage1,lineage2):
        count1=0
        count2=0
        for  lineage in sorted_dict:
                if essential.isequal(lineage1,lineage):
                    count1=sorted_dict[lineage]
                if essential.isequal(lineage2,lineage):
                    count2=sorted_dict[lineage]
        if count1>count2:
            return lineage1
        elif count1<count2:
            return lineage2
        else:
               return np.random.choice([lineage2,lineage1])
            

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
                                       

def tabulate(tracker,sp,gt,test_dic,lineage1,lineage2):
        red= readWrite.readWrite()
        

        #visited_list=[]
        visited,value_list1,value_list,dic=account_trios(tracker)
        

        
    
        dic=dict(sorted(dic.items(),key=lambda item: item[1],reverse=True))
        dic =merge_moving(dic,red)
        
        moving_taxas=querry_lineage(dic,lineage1,lineage2)
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
                                                
                             
                             
              

        return overall_return
        


      

def _parse_list(s: str):
    return [x for x in re.split(r'[,\s/]+', s.strip()) if x]

def parse1():
    parser = argparse.ArgumentParser(description="DAFT Tranform")
    parser.add_argument('--sp', type=str, help="Species tree")
    parser.add_argument(
        '--lineages',
        type=lambda s: s.split('/'),
        help="Test lineage separated  by '/' (e.g. l1/l2)"
    )
    parser.add_argument('--gt_list', nargs='+', help="List of gene trees (each a Newick string)")
    parser.add_argument('--output', type=str, help="Name of output file")
    
    args = parser.parse_args()
    return args


parser = parse1()
sp_string = parser.sp
lineages = parser.lineages 
output =parser.output

lineage1= lineages[0]
lineage2= lineages[1]

output= output+'_'+lineage1[:-1]+'_'+lineage2[:-1]
lis=[]
#sp_string='((((1,2),3),((4,5),6)),7);'
#sp_string='((((((((1,2),3),4),5),6),7),(8,9)));'
gt_list=parser.gt_list



essential= daft_essential()
reco =reconcils()
red= readWrite.readWrite()
write_intro={'idx':[],'Replicate':[],'Path':[],'From_Where_moved':[],'Sibling':[],'What_moved':[],'NNI':[]}
write_intro1={'Moves':[],'Replicate':[],'NNI':[]}



#xit()
#exit()
#exit()
sp_large= red.parse_bio(sp_string)
sp_large.label_internal()
#red.write_introgression(sp_large)
done=[]



for k,gene_tree in enumerate(gt_list):
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
        reco.L_cost=  2
        reco.D_cost=  20000
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
                with ThreadPoolExecutor(max_workers=num_threads) as executor:
                        future = executor.submit(call_reconcILS_function)

                try:
                        li , introgression, tracker,test_dic= future.result(timeout=timeout_duration)
                except TimeoutError:
                        print(f"The function call timed out after {timeout_duration} seconds.")
                        exit()
                        li=[]
                
                
                #print('-----------------------------------------------------Done With reconcils-----------------------------------------------------------------------------')

                if len(tracker)>=1:
                        table=tabulate(tracker,sp,tr,test_dic,lineage1,lineage2)
                        added=0
                        
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
                
                