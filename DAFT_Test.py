import numpy as np
import os, sys
sys.path.append("./DAFT_utils/reconcILS")
sys.path.append("./DAFT_utils/")
from reconcILS import *
from utils_reconcILS import *
import pandas as pd 
import argparse
from DAFT_essential import *
import copy


'''
def render_val(col, v):
    if col in ('Z-value-uncle', 'Z-value-sibling'):
        return '-' if pd.isna(v) else f"{float(v):.2f}"
    elif col in ('total_count', 'uncle_count', 'sibling_count'):
        if pd.isna(v):
            return '-'
        try:
            fv = float(v)
            return str(int(fv)) if fv.is_integer() else str(fv)
        except Exception:
            return str(v)
    else:
        return '' if v is None or (isinstance(v, float) and math.isnan(v)) else str(v)



def _fix_group_nested(df):
    what_list = df['What_moved'].tolist()
    for idx, cu in df['comparison_uncle'].items():
        for wm in what_list:
            if (not isequal(cu, wm)) and isequal_set(cu, wm):
                df.at[idx, 'comparison_uncle'] = wm
                break  
    for idx, cu in df['comparison_sibling'].items():
        for wm in what_list:
            if (not isequal(cu, wm)) and isequal_set(cu, wm):
                df.at[idx, 'comparison_sibling'] = wm
                break  
    return df

def _uncle_count(row):
    
    mask1 = grouped['What_moved'].apply(lambda x: isequal_set(x, row['comparison_uncle'])) 
    mask2 = grouped['Where_at'].apply(lambda x: isequal_set(x, row['Where_at']))
    
    if mask2.any() and not mask1.any():
        return pd.NA          
    return grouped.loc[mask1&mask2, 'total_count'].sum()


def _sibling_count(row):
    
    mask1 = grouped['What_moved'].apply(lambda x: isequal_set(x, row['comparison_sibling'])) 
    mask2 = grouped['Where_at'].apply(lambda x: isequal_set(x, row['Where_at']))
    
    if mask2.any() and not mask1.any():
        return pd.NA          
    return grouped.loc[mask1&mask2, 'total_count'].sum()
'''
def filter_data(data,lineage1,lineage2):
    filtered_data=[]

    red= readWrite.readWrite()
    l1 = red.parse(lineage1)
    l2= red.parse(lineage2)
    l1.label_internal()
    l2.label_internal()
    taxa_1=l1.taxa
    taxa_2=l2.taxa
    for gt in data:
        ##print(gt)
        gt =gt[0]
        ##print(gt)
        gt =gt.replace('e-', '0')
        gt_tr= red.parse(gt)
        gt_tr.label_internal()


        if  essential.current_address(taxa_1,gt_tr) and essential.current_address(taxa_2,gt_tr):
            filtered_data.append([gt])
    
    return filtered_data

def extract_truth(demography):
    text = open(demography).read()

    
    pattern = r"│\s*(\d+)\s*│\s*([A-Z])\s*│"


    matches_id = re.findall(pattern, text)
    id_to_name = {int(match[0]): match[1] for match in matches_id}
    
    mass_migrations = re.findall(
        r"source\s*=\s*(\d+)\s*,\s*dest\s*=\s*(\d+)",
        text
    )
    mass_migrations = {(id_to_name[int(src)], id_to_name[int(dst)]) for src, dst in mass_migrations}
    ##print(mass_migrations)

    return mass_migrations



def write_results(result_sig,out,file="overall_result.csv"):
    df = pd.DataFrame(result_sig)
    #file = "overall_result.csv"

    if os.path.exists(file):
        existing = pd.read_csv(file, index_col=0)
        #df.index = range(existing.index.max() + 1, existing.index.max() + 1 + len(df))
        df.index = [int(out)] * max(len(df),1)

        
        header = False
    else:
        header = True

    df.to_csv(file, mode="a", header=header)


def extract_block(text, start_key):
 
    pattern = rf"{re.escape(start_key)}(.*?)(\*{{5,}})"
    match = re.search(pattern, text, re.S)
    return match.group(1).strip() if match else ""


def parse_inferred_relations(block):
    rows = []

    for line in block.splitlines():
        line = line.strip()
        if not line or "=" in line:
            continue

        # Split before RECEIVER
        if "RECEIVER:" not in line:
            continue

        left, right = line.split("RECEIVER:", 1)
        significant_pairs = left.strip()

        # Extract RECEIVER and DONOR
        receiver = None
        donor = None

        if "AND Donor:" in right:
            receiver_part, donor_part = right.split("AND Donor:", 1)
            receiver = receiver_part.strip()
            donor = donor_part.strip()

        rows.append([
            significant_pairs,
            receiver,
            donor
        ])

    return rows

def extract_direction(out):

    text = open("./DAFT_Direction_"+out+".txt", "r").read()
    ir_block = extract_block(text, "INFERRED RELATIONS")
    df_ir = parse_inferred_relations(ir_block)
    list_return=[]
    for i in df_ir:
        list_return+=[(i[1],i[2])]
        
    return list_return



    

        


def call_direction(sorted_grouped,gene_treefile,sp,out):
    
    list_unique= filter_direction(sorted_grouped,-1.96,5,0,out)
    #print(list_unique)
    #exit()

    sibling_list= essential.find_sibling(sp,[])
   
    list_non_unique=add_sibling_bidirection(sorted_grouped,list_unique,sibling_list)

    
    if len(list_unique)==0:
        print('Zero Significant')

    else:
        #list_non_unique=pairs1

        script = "./DAFT_Direction.py"
        sp = sp_string
        output =out

        argv = [
            sys.executable,         
            script,
            "--sp", str(sp),
            "--gt", str(gene_treefile),
            "--lineages", repr(list_unique),  
            "--lineagesN", repr(list_non_unique),  
            "--output", output,
        ]
        status = os.spawnv(os.P_WAIT, sys.executable, argv)
        exit_code = os.WEXITSTATUS(status) if hasattr(os, "WEXITSTATUS") else (status >> 8)
        if exit_code != 0:
            raise RuntimeError(f"DAFT_Direction.py failed with code {exit_code}")
        
        #exit()
    




def find_sib_lineage_pair(data):
    sib_lineage={}
    for g_t in data:
        #print(g_t)
        g_t[0] = g_t[0].replace('e-', '0')
        tr= red.parse(g_t[0])
        tr.label_internal()
        list_gt_sibling=essential.find_sibling(tr,[])


        for gt_sib in list_gt_sibling:
                            #for lin in newick_lineage:
                            #if isequal(gt_sib[0],lin):
                            #for lin1 in newick_lineage:
                            #if isequal(gt_sib[1],lin1):
                            if (gt_sib[0],gt_sib[1]) in sib_lineage:
                                sib_lineage[(gt_sib[0],gt_sib[1])]+=1
                            elif (gt_sib[1],gt_sib[0]) in sib_lineage:
                                sib_lineage[(gt_sib[1],gt_sib[0])]+=1
                            else:
                                sib_lineage[(gt_sib[0],gt_sib[1])]=1
    return sib_lineage

def accounting(data,sib_lineage):
    result={'Where_at':[],'What_moved':[],'NNI_sp':[],'total_count':[]}
    #result1={'Where_at':[],'What_moved':[],'NNI_sp':[],'total_count':[]}
    visited_dic_set={}
    visited_set_dict ={}
    visited_set_index=[]
    running_index=0
    for ke,valu in sib_lineage.items():
        #print(ke)
        lineage,sibling=ke
        key1_t = red.parse(lineage)
        key_t= red.parse(sibling)
        key1_t.label_internal()
        key_t.label_internal()
        taxa_1=key1_t.taxa
        taxa_2=key_t.taxa
        taxa_1_B,taxa_2_B=taxa_1,taxa_2
        

        if not  essential.current_address(taxa_1,sp) or not essential.current_address(taxa_2,sp):
            #print(lineage,sibling)
            #exit()
            '''
            result1['Where_at']+=[key1_t.to_newick()]
            result1['What_moved']+=[key_t.to_newick()]
            result1['NNI_sp']+=[-1]
            result1['total_count']+=[valu]

            result1['Where_at']+=[key1_t.to_newick()]
            result1['What_moved']+=[key_t.to_newick()]
            result1['NNI_sp']+=[-1]
            result1['total_count']+=[valu]
            pd.DataFrame(result1).to_csv('rev_all.csv',index=False)
            '''
            continue
        if key1_t.isLeaf:
            taxa_1_B={taxa_1}

        if key_t.isLeaf:
            taxa_2_B={taxa_2}


        merged_lineage=frozenset(sorted(taxa_1_B))
        merged_sibling=frozenset(sorted(taxa_2_B))

        if merged_lineage in visited_set_dict and merged_sibling in visited_set_dict:
        #if frozenset(taxa_1_B) in visited_set and frozenset(taxa_2_B) in visited_set:
            key_string_lineage=visited_set_dict[merged_lineage]
            key_string_sibling=visited_set_dict[merged_sibling]
            if sorted([key_string_lineage,key_string_sibling]) in visited_set_index:
                #dist=findDist(sp,key_t.taxa,key1_t.taxa)-2
                #if dist>=0:
                #print(sorted([key_string_lineage,key_string_sibling]), visited_set_index)
                #print(visited_set_dict)
                #print(sorted([key_string_lineage,key_string_sibling]) in visited_set_index)
                #print(lineage,sibling,taxa_1,taxa_2)
                #index_lineage = {iii for iii, iere in enumerate(result['Where_at']) if iere == lineage}
                #index_sibling = {iii for iii, iere in enumerate(result['What_moved']) if iere == sibling}

                index_lineage=set()
                for iii, iere in enumerate(result['Where_at']):
                        iere_t= red.parse(iere)
                        iere_t.label_internal()
                        #print(iere_t,taxa_1)
                        if iere_t.taxa==taxa_1:
                            index_lineage.add(iii)

                index_sibling=set()
                for iii, iere in enumerate(result['What_moved']):
                        iere_t= red.parse(iere)
                        iere_t.label_internal()
                        if iere_t.taxa==taxa_2:
                            index_sibling.add(iii)

                #print(result)
                #print(''.join(sorted(''.join(taxa_1_B))),''.join(sorted(''.join(taxa_2_B))),index_lineage,index_sibling)
                target_index = index_lineage.intersection(index_sibling).pop()
                #print(list(target_index)[0])

                result['total_count'][target_index]+=valu

                index_lineage=set()
                for iii, iere in enumerate(result['What_moved']):
                        iere_t= red.parse(iere)
                        iere_t.label_internal()
                        if iere_t.taxa==taxa_1:
                            index_lineage.add(iii)

                index_sibling=set()
                for iii, iere in enumerate(result['Where_at']):
                        iere_t= red.parse(iere)
                        iere_t.label_internal()
                        if iere_t.taxa==taxa_2:
                            index_sibling.add(iii)
                target_index_reverse = index_lineage.intersection(index_sibling).pop()
                
                result['total_count'][target_index_reverse]+=valu
                pd.DataFrame(result).to_csv('rev_n.csv',index=False)

            else:
                dist=essential.findDist(sp,key_t.taxa,key1_t.taxa)-2
                
                if dist>=0:
                    lineage= visited_dic_set[visited_set_dict[merged_lineage]]
                    sibling= visited_dic_set[visited_set_dict[merged_sibling]]
                    
                    #visited_set.add(frozenset(taxa_1))
                    #visited_set.add(frozenset(taxa_2))
                    
                    #visited_set_dict[''.join(sorted(''.join(taxa_1)))]=running_index
                    #visited_set_dict[''.join(sorted(''.join(taxa_2)))]=running_index+1
                    visited_set_index+=[sorted([key_string_lineage,key_string_sibling])]
                    result['Where_at']+=[lineage]
                    result['What_moved']+=[sibling]
                    result['NNI_sp']+=[dist]
                    result['total_count']+=[valu]

                    result['Where_at']+=[sibling]
                    result['What_moved']+=[lineage]
                    result['NNI_sp']+=[dist]
                    result['total_count']+=[valu]
                    #print('flag1')
                    #running_index+=2   
                    pd.DataFrame(result).to_csv('rev_n.csv',index=False)

        else:
            dist=essential.findDist(sp,key_t.taxa,key1_t.taxa)-2
            if dist>=0:
                #print('flag2',running_index)

                if merged_lineage not in visited_set_dict:
                    visited_set_dict[merged_lineage]=running_index
                    visited_dic_set[running_index]=lineage
                    running_index+=1
                else:
                    lineage= visited_dic_set[visited_set_dict[merged_lineage]]
                
                if  merged_sibling not in visited_set_dict:
                    visited_set_dict[merged_sibling]=running_index
                    visited_dic_set[running_index]=sibling
                    running_index+=1
                else:
                    sibling= visited_dic_set[visited_set_dict[merged_sibling]]


                visited_set_index+=[sorted([ visited_set_dict[merged_lineage],visited_set_dict[merged_sibling]])]
                
                result['Where_at']+=[lineage]
                result['What_moved']+=[sibling]
                result['NNI_sp']+=[dist]
                result['total_count']+=[valu]

                result['Where_at']+=[sibling]
                result['What_moved']+=[lineage]
                result['NNI_sp']+=[dist]
                result['total_count']+=[valu]

                #print('flag2',running_index)


                pd.DataFrame(result).to_csv('rev_n.csv',index=False)  
    return result 

def sorting_arrangement():
    df = pd.read_csv("./rev_n.csv", sep=',')
    grouped = df.groupby(['Where_at', 'What_moved', 'NNI_sp'], as_index=False)['total_count'].sum()
    sorted_grouped = grouped.sort_values(by=['Where_at', 'NNI_sp'])
    sp_tree_lineages= essential.find_all_lineage(sp)
                
    #print(sorted_grouped)

    sorted_grouped["Where_at"] = sorted_grouped["Where_at"].apply(
        lambda v: essential.match_lineage(v, sp_tree_lineages)
    )
    sorted_grouped["What_moved"] = sorted_grouped["What_moved"].apply(
        lambda v: essential.match_lineage(v, sp_tree_lineages)
    )
    return grouped,sorted_grouped,df

def address_to_newick(list_sp_tree):
    newick_lineage=[]
    for i in list_sp_tree:
        newick_lineage.append(i.to_newick())
    return newick_lineage

def add_sibling_bidirection(sorted_grouped,pairs,sibling_list):
    pairs1=[]
    for idx,na in enumerate(pairs):
        wm,wa=na[0],na[1]
        for sibi in sibling_list:
            sibi1= sibi[0]
            sibi2= sibi[1]
            if essential.isequal_set(wm,sibi1) and not essential.is_in(wa,sibi2):
                pairs1.append(([wm,sibi2],wa))

            if essential.isequal_set(wm,sibi2) and not essential.is_in(wa,sibi1):
                pairs1.append(([wm,sibi1],wa))
            
            if essential.isequal_set(wa,sibi1) and not essential.is_in(wm,sibi2):
                pairs1.append(([wa,sibi2],wm))

            if essential.isequal_set(wa,sibi2) and not essential.is_in(wm,sibi1):
                pairs1.append(([wa,sibi1],wm))
    return pairs1

def switch_uncle(lineage,uncle_list,where_at):
    #print(uncle_list)
    for pair in uncle_list:

        key=pair[0]
        uncle =pair[1]
        #print(key,uncle)
        #print(lineage,uncle,where_at,essential.is_in(where_at,uncle))
    
        if essential.isequal_set(lineage,uncle) and not essential.is_in(where_at,key):
            return key
            
    return None

def find_distance_sib(What_moved,comp_uncle,where_at):
    if comp_uncle:
        Tree1= red.parse(What_moved)
        Tree2= red.parse(comp_uncle)
        Tree3= red.parse(where_at)

        Tree1.label_internal()
        Tree2.label_internal()
        Tree3.label_internal()

        distance1=abs(essential.findDist(sp,Tree1.taxa,Tree3.taxa)-2)
        distance2=abs(essential.findDist(sp,Tree2.taxa,Tree3.taxa)-2)
        return distance1,distance2
    else:
        return 0,0


def put_sibling_uncle(grouped,sorted_grouped,df):
    sorted_grouped['comparison_sibling'] = sorted_grouped['What_moved'].apply(lambda x: essential.map_pair_sibling(x,sibling_list))
    #print('--'*100)
    sorted_grouped['comparison_uncle'] = sorted_grouped['What_moved'].apply(lambda x: essential.map_pair_uncle(x,uncle_list))

    #print(uncle_list)
    #for i in zip(sorted_grouped['What_moved'],sorted_grouped['comparison_uncle']):
    #print(i)
    #exit()

    for where_at, group in sorted_grouped.groupby('Where_at',as_index=False):
        for idx, row in group.iterrows():
            What_moved=row['What_moved']
            comp_uncle= row['comparison_uncle']
            comp_sibling=row['comparison_sibling']
            NNI_sp =row['NNI_sp']
            #print(where_at,comp_uncle,What_moved)
            #print('---------------------------------')
            if NNI_sp==0:
                sorted_grouped.loc[idx, 'comparison_sibling'] = None
                sorted_grouped.loc[idx, 'comparison_uncle'] = None
                continue
            
            if essential.isequal_set(comp_sibling,where_at) or essential.is_in(where_at,comp_sibling):
                #print('222222222222222222222222')
                sorted_grouped.loc[idx, 'comparison_sibling'] = None

            if essential.isequal_set(comp_uncle,where_at):
                #print('2222222222222222222222222121212')
                sorted_grouped.loc[idx, 'comparison_uncle'] = None
                continue
            if essential.is_in(where_at,comp_uncle):
                new_uncle = essential.map_pair_sibling(comp_uncle,sibling_list)
                #print(where_at,What_moved, new_uncle,comp_uncle)
                #exit()
                if essential.isequal_set(new_uncle,where_at) or essential.is_in(What_moved,new_uncle):
                    sorted_grouped.loc[idx, 'comparison_uncle'] = None
                    #print(where_at,What_moved, new_uncle,comp_uncle)
               
                    #print('here')
                    #exit()
                    continue
                else:
                    #print('here iss')
                    comp_uncle =new_uncle
                    sorted_grouped.loc[idx, 'comparison_uncle']=comp_uncle
                    continue
                
            ##print(What_moved,comp_uncle)
            if comp_uncle:
                distance1, distance2= find_distance_sib(What_moved,comp_uncle,where_at)
            

                ##print('None',What_moved,comp_uncle,where_at,distance1,distance2)
                if distance1==0 or distance2==0:
                    sorted_grouped.loc[idx, 'comparison_uncle'] = None
                elif  distance2>distance1:
                    new_uncle = switch_uncle(What_moved,uncle_list,where_at)
                    #print(new_uncle)
                    if essential.isequal_set(new_uncle,where_at) or essential.is_in(where_at,new_uncle):
                        sorted_grouped.loc[idx, 'comparison_uncle'] = None
                    else:
                        distance1, distance2= find_distance_sib(What_moved,new_uncle,where_at)
                        if distance1==0 or distance2==0:
                            sorted_grouped.loc[idx, 'comparison_uncle'] = None
                        elif distance2>distance1:
                            sorted_grouped.loc[idx, 'comparison_uncle'] = None
                        else:
                            sorted_grouped.loc[idx, 'comparison_uncle'] = new_uncle
                    continue
                    

    return sorted_grouped

def tuple_equal(t1, t2):
    a, b = t1
    c, d = t2
    
    a_tree =red.parse(a)
    b_tree =red.parse(b)
    c_tree =red.parse(c)
    d_tree =red.parse(d)
    
    a_tree.label_internal()
    b_tree.label_internal()
    c_tree.label_internal()
    d_tree.label_internal()
    
    #a=a_tree.taxa
    #b=b_tree.taxa
    #c=c_tree.taxa
    #d=d_tree.taxa
    ##print(a,b,c,d)
    
    
    return (essential.isequal(a, c) and essential.isequal(b, d)) or \
           (essential.isequal(a, d) and essential.isequal(b, c))



def convert_species(df):
    species_to_letters = {
    'Carlitosyrichta': 'A',
    'Cebuscapucinus': 'B',
    'Cercocebusatys': 'C',
    'Macacafascicularis': 'D',
    'Macacanemestrina': 'E',
    'Theropithecusgelada': 'F',
    'Papioanubis': 'G',
    'Macacamulatta': 'H',
    'Mandrillusleucophaeus': 'I',
    'Chlorocebussabaeus': 'J',
    'Colobusangolensis': 'K',
    'Piliocolobustephrosceles': 'L',
    'Rhinopithecusbieti': 'M',
    'Rhinopithecusroxellana': 'N',
    'Gorillagorilla': 'O',
    'Homosapiens': 'P',
    'Panpaniscus': 'Q',
    'Pantroglodytes': 'R',
    'Pongoabelii': 'S',
    'Nomascusleucogenys': 'T',
    'Saimiriboliviensis': 'U',
    'Aotusnancymaae': 'V',
    'Callithrixjacchus': 'W',
    }
    #letter_to_species ={value:key for key,value in species_to_letters.keys()}
    
    for idx,row in df.iterrows():
        old_what_moved = row['Test_lineage']
        old_where_at = row['Focal_lineage']
        old_comp_sib = row['comparison_sibling']
        old_comp_uncle = row['comparison_uncle']
        
        old_what_moved_tree = red.parse(old_what_moved)
        old_where_at_tree = red.parse(old_where_at)
        
        new_what_moved =old_what_moved_tree.to_newick_change(species_to_letters)
        new_where_at = old_where_at_tree.to_newick_change(species_to_letters)
        
        new_comp_sib=old_comp_sib
        new_comp_uncle= old_comp_uncle
        if old_comp_sib:
            old_comp_sib_tree =red.parse(old_comp_sib)
            
            new_comp_sib = old_comp_sib_tree.to_newick_change(species_to_letters)
            
        if old_comp_uncle:
            old_comp_uncle_tree =red.parse(old_comp_uncle)
            new_comp_uncle = old_comp_uncle_tree.to_newick_change(species_to_letters)
        

        
        
        
        df.loc[idx,'Test_lineage']=new_what_moved
        df.loc[idx,'Focal_lineage']=new_where_at
        df.loc[idx,'comparison_sibling']=new_comp_sib
        df.loc[idx,'comparison_uncle']=new_comp_uncle
    return df
        
    
def find_level_from_below(address):
    hop=0
    
    stack=[(address,hop)]
    #li=[]
    
    while stack:
        current_node,hop1 =stack.pop()
        hop =max(hop1,hop)
        if current_node.isLeaf:
            continue
        stack.append((current_node.leftChild,hop+1))
        stack.append((current_node.rightChild,hop+1))
    return hop-1
        
        
def find_level_from_above(address):
    hop=0
    
    stack=[(address,hop)]
    #li=[]

    
    while stack:
        current_node,hop =stack.pop()
        #hop =max(hop1,hop)
        if current_node:
            if current_node.parent:
                stack.append((current_node.parent,hop+1))
        
    return hop+1
        
    

def filter_direction(sorted_grouped,Z,tc,nni,out):
    filtered_df_1 = sorted_grouped[
        (sorted_grouped['Z-value-uncle'] <= Z)&(sorted_grouped['NNI_sp']>nni)&((sorted_grouped['total_count']+sorted_grouped['uncle_count'])>tc)
    ].assign(flag=1)
    filtered_df_2 = sorted_grouped[
        (sorted_grouped['Z-value-sibling'] <= Z)&(sorted_grouped['NNI_sp']>nni)&((sorted_grouped['total_count']+sorted_grouped['sibling_count'])>tc)
    ].assign(flag=2)
    
    ##print('hellp')
    
    concat_df = pd.concat([filtered_df_1, filtered_df_2], ignore_index=True)
    ##print('hellp222')
    
    filtered_df_3=concat_df.copy()
    filtered_df_3.rename(columns={'What_moved': 'Test_lineage', 'Where_at': 'Focal_lineage','total_count':'Test_count'}, inplace=True)
    #filtered_df_3.to_csv('important_'+out+'.csv',index=False)
    #write_results(filtered_df_3.to_dict(),file="Important.csv")
    
    filtered_df_3["idx"] = out 
    #Focal_lineage,Test_lineage,NNI_sp,Test_count,comparison_sibling,comparison_uncle,sibling_count,Z-value-sibling,uncle_count,Z-value-uncle,flag

    expected_columns = [
    "Focal_lineage",
    "Test_lineage",
    "NNI_sp",
    "Test_count",
    "comparison_sibling",
    "comparison_uncle",
    "sibling_count",
    "Z-value-sibling",
    "uncle_count",
    "Z-value-uncle",
    "flag",
    "idx"
    ]

    #filtered_df_3 = filtered_df_3.copy()
    filtered_df_3["idx"] = out

    if filtered_df_3.empty:
        filtered_df_3 = pd.DataFrame(columns=expected_columns)
        filtered_df_3.loc[0, "idx"] = out


    
    #filtered_df_3= convert_species(filtered_df_3)
    if os.path.exists("./Summary.csv"):
        filtered_df_3.to_csv("./Summary.csv", mode="a", header=False, index=False)
    else:
        filtered_df_3.to_csv("./Summary.csv", mode="w", header=True, index=False)
        
        
    pairs = list(zip(concat_df['What_moved'], concat_df['Where_at']))
    ##print(pairs)
    pairs1= []
    not_include=[]
    for p1 in range(len(pairs)):
        if p1 not in not_include:
            pairs1+=[pairs[p1]]
        for p2 in range(p1,len(pairs)):
            if tuple_equal(pairs[p1],pairs[p2]):
                not_include.append(p2)
                continue
                
    ##print(pairs1)
    #exit()
    return pairs1


def put_z(sorted_grouped,grouped):
    
    sibling_counts = (
        sorted_grouped.apply(essential._sibling_count, axis=1,args=(grouped,))
    )
    
    uncle_counts = (
        sorted_grouped.apply(essential._uncle_count, axis=1,args=(grouped,))
    )
    
    sibling_counts = sibling_counts.fillna(0) 
    sibling_counts = sibling_counts.astype(float)
    
    uncle_counts = uncle_counts.fillna(0) 
    uncle_counts = uncle_counts.astype(float)
    
    sorted_grouped['sibling_count'] = sibling_counts
    sorted_grouped['uncle_count'] = uncle_counts
    

    
    essential.compute_z(sorted_grouped,'sibling',df)
    
    


    essential.compute_z(sorted_grouped,'uncle',df)
    sorted_grouped = (
        sorted_grouped
        .groupby('Where_at', group_keys=False, dropna=False)
        .apply(essential._fix_group_nested)
    )

    return sorted_grouped


def write_significance(sorted_grouped,out_file,labeled_sp):
    sorted_grouped = sorted_grouped.sort_values(by=[ 'NNI_sp'])

    cols = [
        'NNI_sp',
        'Test_lineage',
        'Test_count',
        'comparison_uncle',
        'uncle_count',
        'Z-value-uncle',
        'comparison_sibling',
        'sibling_count',
        'Z-value-sibling',
    ]

    if sibling_flag == '0':
        cols = [
            'NNI_sp',
            'Test_lineage',
            'Test_count',
            'comparison_uncle',
            'uncle_count',
            'Z-value-uncle'
        ]

    sorted_grouped.rename(columns={"What_moved": "Test_lineage"}, inplace=True)
    sorted_grouped.rename(columns={"total_count": "Test_count"}, inplace=True)
    sorted_grouped.rename(columns={"Where_at": "Focal_lineage"}, inplace=True)
    #sorted_grouped =convert_species(sorted_grouped)

    padding = 2 
    widths = {}
    for col in cols:
        header_len = len(col)
        if col in sorted_grouped.columns:
            data_lens = sorted_grouped[col].map(lambda v: len(essential.render_val(col, v)))
        else:
            data_lens = pd.Series([0])
        widths[col] = max([header_len] + data_lens.tolist())

    # --- separators ---
    SPECIAL_AFTER = ['Test_count', 'Z-value-uncle']  # put 3 tabs after these columns
    TAB_STR = '\t' * 5
    SPACE_SEP = ' ' * padding
    TABSTOP = 8  

    with open('./DAFT_Significance_'+out_file+'.txt', 'w') as oe:
        oe.write(f"Species Tree = {essential.to_newick_with_id(labeled_sp)}\n")
        oe.write("=" * 40 + "\n")

        for where_at, group in sorted_grouped.groupby('Focal_lineage', as_index=False):
            oe.write(f"Focal_lineage = {where_at}\n")

            
            header = ''
            for i, col in enumerate(cols, 1):  
                header += f"{col:<{widths[col]}}"
                if i < len(cols):  
                    header += TAB_STR if col in SPECIAL_AFTER else SPACE_SEP
            header += "\n"

            oe.write(header)
            oe.write("=" * len(header.expandtabs(TABSTOP).rstrip()) + "\n")

            
            for _, row in group.iterrows():
                line = ''
                for i, col in enumerate(cols, 1):
                    val = row[col] if col in row else np.nan
                    s = essential.render_val(col, val)
                    line += f"{s:<{widths[col]}}"
                    if i < len(cols):
                        
                        line += TAB_STR if col in SPECIAL_AFTER else SPACE_SEP
                line += "\n"
                oe.write(line)

            oe.write("-" * len(header.expandtabs(TABSTOP).rstrip()) + "\n")


def convert_excel(out):
        script = "./DAFT_produce_excel.py"
        argv = [
                sys.executable,         
                script,
                "--output",out,
                
            ]
        status = os.spawnv(os.P_WAIT, sys.executable, argv)
        exit_code = os.WEXITSTATUS(status) if hasattr(os, "WEXITSTATUS") else (status >> 8)
        if exit_code != 0:
            raise RuntimeError(f"DAFT_produce_excel.py failed with code {exit_code}")

def clean_folder(out):
        script = "./clean_folder.py"
        argv = [
                sys.executable,         
                script,
                "--output",out,
            ]
        status = os.spawnv(os.P_WAIT, sys.executable, argv)
        exit_code = os.WEXITSTATUS(status) if hasattr(os, "WEXITSTATUS") else (status >> 8)
        if exit_code != 0:
            raise RuntimeError(f"clean_folder.py failed with code {exit_code}")

def is_leaf(x):
    if x:
        #print(x)
        T1= red.parse(x)
        if T1.isLeaf:
            return True
        else:
            return False
    return False
    
def extract_leaves(x):
    if x:
        T1= red.parse(x)
        return [T1.leftChild.to_newick(),T1.rightChild.to_newick()]
    return None

def correct_value(x,leaf_grouped,where_at):
    for idx, row in leaf_grouped.iterrows():
        print(row)
    
    pass


def extract_mrca_linage(sp_string,lin1,focal):
    sp= red.parse(sp_string)
    sp.label_internal()
    old_1 =red.parse(lin1)
    old_1.label_internal()
    
    old_2 =red.parse(focal)
    old_2.label_internal()
    curr_1=  essential.current_address(old_1.taxa,sp)
    
    
    curr_2=  essential.current_address(old_2.taxa,sp)
    curr_2_taxa =old_2.taxa
    if curr_2.isLeaf:
        curr_2_taxa= {old_2.taxa}
    stack=[curr_1]
    li=[]
    hop=0
    
    while stack:
        current_node =stack.pop()
        if current_node:
            if curr_2_taxa.issubset(current_node.taxa):
                return li 
            else:
                stack.append(current_node.parent)
                hop+=1
                if hop>0:
                    print(old_1.taxa,old_2.taxa,current_node.to_newick(),current_node.taxa)
                    sib_1= [aer  for aer in current_node.parent.children if aer !=current_node][0]
                    #if sib_1.taxa != old_2.taxa:
                    li.append(sib_1.to_newick())
                    if not sib_1.isLeaf:
                        
                        leaves= extract_leaves(sib_1.to_newick())
                        li+=leaves

                    #li+=[lim for lim in leaves if lim if not essential.isequal(lim,focal)]
    return li
        

def correct_accounting_by_product(data,leaf1,where_at,leaf2,uncle,sp_string):
    
    new_gt_2 = '(('+leaf1+','+where_at+'),'+uncle+')'
    new_gt_2_tree = red.parse(new_gt_2)
    new_gt_2_tree.label_internal()
    
    new_gt_3 = '(('+leaf2+','+where_at+'),'+uncle+')'
    new_gt_3_tree = red.parse(new_gt_3)
    new_gt_3_tree.label_internal()
    #print(new_gt_1)
    
    print('/;/;;;;;;',leaf1,where_at,leaf2,uncle, new_gt_2,new_gt_3)
    
    #topo_X='(('+leaf1+','+where_at+'),'+leaf2+');'
    #topo_y='(('+leaf2+','+where_at+'),'+leaf1+');'
    #topo_z='(('+leaf2+','+leaf1+'),'+where_at+');'

    count_k=0
    count_l=0
    
    count_k_1=0
    
    count_l_1=0
    #accout= {'gt':[],'topo':[]}
    sp = red.parse(sp_string)
    sp.label_internal()
    #print(sp,where_at,uncle)
    dist= essential.find_dist_string(sp,where_at,uncle)
    if dist <0:
        return count_l,count_k
    
    L2=red.parse(leaf2)
    L1=red.parse(leaf1)
    for g_t in data:
        #print(g_t)
        g_t[0] = g_t[0].replace('e-', '0')
        tr= red.parse(g_t[0])
        tr.label_internal()
        
        
    
        addr2= essential.current_address(new_gt_2_tree.taxa,tr)
        addr3= essential.current_address(new_gt_3_tree.taxa,tr)
        
        #print(addr1,tr.to_newick(),new_gt_1_tree.taxa)
        #print(topo_X)
        #print(where_at,'('+leaf1+','+where_at+')')
        
        
        if addr2:
            count_k_1+=1
            
            #print('add2',addr2.to_newick(),leaf1,where_at,leaf2,uncle)
            if (essential.isequal_set(addr2.leftChild.to_newick(),'('+leaf1+','+where_at+')')  and  essential.isequal_set(addr2.rightChild.to_newick(),uncle)) or (essential.isequal_set(addr2.rightChild.to_newick(),'('+leaf1+','+where_at+')')  and  essential.isequal_set(addr2.leftChild.to_newick(),uncle)):  
                print(addr2.to_newick())
                count_k+=1
                #accout['gt']+=[tr.to_newick()]
                #accout['topo']+=[addr2.to_newick()]
                
            
        
        
        if addr3:
            count_l_1+=1
            
            #print('add3',addr3.to_newick(),leaf1,where_at,leaf2,uncle)
            if (essential.isequal_set(addr3.leftChild.to_newick(),'('+leaf2+','+where_at+')')  and  essential.isequal_set(addr3.rightChild.to_newick(),uncle)) or (essential.isequal_set(addr3.rightChild.to_newick(),'('+leaf2+','+where_at+')')  and  essential.isequal_set(addr3.leftChild.to_newick(),uncle)):  
                print(addr3.to_newick())
                count_l+=1
                #accout['gt']+=[tr.to_newick()]
                #accout['topo']+=[addr3.to_newick()]
                



                
            
    
        
    return count_l,count_k#,
'''


def correct_accounting_by_product(data,leaf1,where_at,leaf2,uncle,sp_string):
    
    new_gt_2 = '((('+leaf1+','+where_at+'),'+uncle+'),'+leaf2+')'
    new_gt_2_tree = red.parse(new_gt_2)
    new_gt_2_tree.label_internal()
    
    new_gt_3 = '(('+leaf2+','+where_at+'),'+uncle+')'
    new_gt_3_tree = red.parse(new_gt_3)
    new_gt_3_tree.label_internal()
    #print(new_gt_1)
    
    print(new_gt_2,new_gt_3)
    
    #topo_X='(('+leaf1+','+where_at+'),'+leaf2+');'
    #topo_y='(('+leaf2+','+where_at+'),'+leaf1+');'
    #topo_z='(('+leaf2+','+leaf1+'),'+where_at+');'

    count_k=0
    count_l=0
    
    count_k_1=0
    
    count_l_1=0
    #accout= {'gt':[],'topo':[]}
    sp = red.parse(sp_string)
    sp.label_internal()
    #print(sp,where_at,uncle)
    dist= essential.find_dist_string(sp,where_at,uncle)
    if dist <0:
        return count_l,count_k
    
    L2=red.parse(leaf2)
    L1=red.parse(leaf1)
    for g_t in data:
        #print(g_t)
        g_t[0] = g_t[0].replace('e-', '0')
        tr= red.parse(g_t[0])
        tr.label_internal()
        
        
    
        addr2= essential.current_address(new_gt_2_tree.taxa,tr)
        #addr3= essential.current_address(new_gt_3_tree.taxa,tr)
        
        #print(addr1,tr.to_newick(),new_gt_1_tree.taxa)
        #print(topo_X)
        #print(where_at,'('+leaf1+','+where_at+')')
        
        
        if addr2:
            #count_k_1+=1
            
            #print('add2',addr2.to_newick(),leaf1,where_at,leaf2,uncle)
            if (essential.isequal_set(addr2.leftChild.to_newick(),'('+leaf1+','+where_at+')')  and  essential.isequal_set(addr2.rightChild.to_newick(),'('+leaf2+','+uncle+')')) or (essential.isequal_set(addr2.rightChild.to_newick(),'('+leaf1+','+where_at+')')  and  essential.isequal_set(addr2.leftChild.to_newick(),'('+leaf2+','+uncle+')')):  
                print(addr2.to_newick())
                count_k+=1
                #accout['gt']+=[tr.to_newick()]
                #accout['topo']+=[addr2.to_newick()]
                
            
        
        
            #if addr3:
            #count_l_1+=1
            
            #print('add3',addr3.to_newick(),leaf1,where_at,leaf2,uncle)
            if (essential.isequal_set(addr2.leftChild.to_newick(),'('+leaf2+','+where_at+')')  and  essential.isequal_set(addr2.rightChild.to_newick(),'('+leaf1+','+uncle+')')) or (essential.isequal_set(addr2.rightChild.to_newick(),'('+leaf2+','+where_at+')')  and  essential.isequal_set(addr2.leftChild.to_newick(),'('+leaf1+','+uncle+')')):  
                print(addr2.to_newick())
                count_l+=1
                #accout['gt']+=[tr.to_newick()]
                #accout['topo']+=[addr3.to_newick()]
                



                
            
    
        
    return count_l,count_k#,
''' 

def correct_accounting(data,leaf1,where_at,leaf2,uncle):
    print(leaf1,where_at,leaf2,uncle)
    
    new_gt_1 = '(('+leaf1+','+where_at+'),'+leaf2+')'
    new_gt_1_tree = red.parse(new_gt_1)
    new_gt_1_tree.label_internal()

    new_gt_2 = '(('+leaf1+','+where_at+'),'+uncle+')'
    new_gt_2_tree = red.parse(new_gt_2)
    new_gt_2_tree.label_internal()
    
    new_gt_3 = '(('+leaf2+','+where_at+'),'+uncle+')'
    new_gt_3_tree = red.parse(new_gt_3)
    new_gt_3_tree.label_internal()
    #print(new_gt_1)
    
    #topo_X='(('+leaf1+','+where_at+'),'+leaf2+');'
    #topo_y='(('+leaf2+','+where_at+'),'+leaf1+');'
    #topo_z='(('+leaf2+','+leaf1+'),'+where_at+');'
    
    count_x=0
    count_y=0
    count_z=0
    count_m=0
    print(new_gt_1,new_gt_2)

    accout= {'gt':[],'topo':[]}
        
    for g_t in data:
        #print(g_t)
        g_t[0] = g_t[0].replace('e-', '0')
        tr= red.parse(g_t[0])
        tr.label_internal()
        
        
        addr1= essential.current_address(new_gt_1_tree.taxa,tr)
        
        #print(addr1,tr.to_newick(),new_gt_1_tree.taxa)
        #print(topo_X)
        #print(where_at,'('+leaf1+','+where_at+')')
        if addr1:
            
            count_m+=1
            print(addr1.to_newick())
            #print('add1',addr1.to_newick(),leaf1,where_at,leaf2,uncle)
            if (essential.isequal_set(addr1.leftChild.to_newick(),'('+leaf1+','+where_at+')')  and  essential.isequal_set(addr1.rightChild.to_newick(),leaf2)) or (essential.isequal_set(addr1.rightChild.to_newick(),'('+leaf1+','+where_at+')')  and  essential.isequal_set(addr1.leftChild.to_newick(),leaf2)):  
                #print(addr1.to_newick())
                accout['gt']+=[tr.to_newick()]
                accout['topo']+=[addr1.to_newick()]
                count_x+=1
            elif  (essential.isequal_set(addr1.leftChild.to_newick(),'('+leaf2+','+where_at+')')  and  essential.isequal_set(addr1.rightChild.to_newick(),leaf1)) or (essential.isequal_set(addr1.rightChild.to_newick(),'('+leaf2+','+where_at+')')  and  essential.isequal_set(addr1.leftChild.to_newick(),leaf1)):
                #print(addr1.to_newick())
                accout['gt']+=[tr.to_newick()]
                accout['topo']+=[addr1.to_newick()]
                count_y+=1
            elif (essential.isequal_set(addr1.leftChild.to_newick(),'('+leaf2+','+leaf1+')')  and  essential.isequal_set(addr1.rightChild.to_newick(),where_at)) or (essential.isequal_set(addr1.rightChild.to_newick(),'('+leaf2+','+leaf1+')')  and  essential.isequal_set(addr1.leftChild.to_newick(),where_at)):
                #print('topo_z')
                #accout['gt']+=[tr.to_newick()]
                #accout['topo']+=[addr1.to_newick()]
                count_z+=1
            else:
                efg=1
                #print('you are wrong')
                #print(addr1.to_newick(),new_gt_1,leaf1,where_at,leaf2,uncle)
                #exit()
        
        

    #print(len(accout['gt']))
    #import pprint
    #pprint.pprint(accout)
        
    pd.DataFrame(accout).to_csv('some_accounting.csv',index=False)
            
                
            
    
    count_x= count_m-count_z
    count_y =0
    return count_x,count_y,count_z,new_gt_1_tree.taxa #,count_l,count_k



def extract_all_lineages(tree):
    T1= red.parse(tree)
    
    stack =[T1]
    li=[]
    while stack:
        current_node = stack.pop()
        if current_node.isLeaf:
            li.append(current_node.to_newick())
        
        else:
            li.append(current_node.to_newick())
            stack.append(current_node.leftChild)
            stack.append(current_node.rightChild)
    return li

def remining_lineage(leaf, tree):
    T1 = red.parse(tree)
    T2 = red.parse(leaf)

    target_taxa = T2.taxa
    stack = [T1]

    while stack:
        node = stack.pop()
        if node is None:
            continue

        if node.taxa == target_taxa:
            parent = node.parent
            if parent is None:
                continue

            # find sibling
            sibling = parent.leftChild if parent.rightChild is node else parent.rightChild
            grandparent = parent.parent

            # case: parent is root
            if grandparent is None:
                sibling.parent = None
                return sibling

            # reconnect grandparent -> sibling
            if grandparent.leftChild is parent:
                grandparent.leftChild = sibling
            elif grandparent.rightChild is parent:
                grandparent.rightChild = sibling

            # also keep children list in sync
            grandparent.children = [
                grandparent.leftChild,
                grandparent.rightChild
            ]

            sibling.parent = grandparent

            return T1

        stack.append(node.leftChild)
        stack.append(node.rightChild)

    return T1
    
    
    
    

def recurssive_correction(data,old_clade,where_at,scale,lineage_middle,sp_string,NNI_SP_test,NNI_SP_sib,NNI_SP_uncle,test_count,sibling_test_count,uncle_test_count,flag_4):
    
    leaf1=old_clade.leftChild.to_newick()
    leaf2=old_clade.rightChild.to_newick()
    new_count_sibling=0
    new_count_test=0
    new_count_uncle=0
    
    test_comparision_sibling,test_comparison_uncle=test_count[0],test_count[1]
    if sibling_test_count:
        sibling_comparision_sibling=sibling_test_count
        
    if uncle_test_count:
        uncle_comparision_sibling=uncle_test_count
    '''
    for line in lineage_middle:
        count_l,count_k= correct_accounting_by_product(data,leaf1[:-1],where_at[:-1],leaf2[:-1],line[:-1],sp_string)
        new_count+=scale* (count_l+count_k)
        
        print(count_l,leaf1[:-1],leaf2[:-1] ,count_k,('-------------------------------'))
    '''
    hop=0
    stack =[(old_clade,hop)]
    old_clade_nw=old_clade.to_newick()
    sp = red.parse(sp_string)
    sp.label_internal()
    #print(sp,where_at,uncle)
    
    
                            
                    
    lineage_middle_mrca= extract_mrca_linage(sp_string,old_clade.to_newick()[:-1],where_at[:-1])
    lineage_middle_full= extract_mrca_linage(sp_string,old_clade.to_newick()[:-1],sp_string)

                    
    while stack:
        current_node,hop =stack.pop()
        
        if current_node.isLeaf:
            continue
        else:
            if hop<3:
                stack.append((current_node.leftChild,hop+1))
                stack.append((current_node.rightChild,hop+1))
                
                leaf1= current_node.leftChild.to_newick()
                leaf2= current_node.rightChild.to_newick()

                #if  not essential.isequal(leaf1,old_clade) and not essential.isequal(lf,leaves[0]) and not essential.isequal(lf,leaves[1]):
                
                dist= essential.find_dist_string(sp,current_node.to_newick(),where_at)
                taxa_size =len(current_node.taxa)
                
                if dist>0 or taxa_size>0:
                    if hop!=0:
                        
                        remining_lineage1 = remining_lineage(leaf1,old_clade_nw)
                        remining_lineage2 = remining_lineage(leaf2,old_clade_nw)
                        print(leaf1)
                        print(leaf1,'--232435'*3,remining_lineage1.to_newick())
                        #dist= essential.find_dist_string(sp,current_node.to_newick(),where_at)
                        #for line in lineage_middle:
                        count_l,count_k,kiase,taxa_mrca= correct_accounting(data,leaf1[:-1],where_at[:-1],remining_lineage1.to_newick()[:-1],sp_string)
                    
                        print(leaf1[:-1],leaf2[:-1] ,count_l,count_k,('-------------------------------'))
                        new_count_sibling +=  1.5* (count_l+count_k)
                        new_count_test +=  1.5* (count_l+count_k)
                        new_count_uncle +=  1.5* (count_l+count_k)
                        

               
                        count_l,count_k,kiase,taxa_mrca= correct_accounting(data,leaf2[:-1],where_at[:-1],remining_lineage2.to_newick()[:-1],sp_string)
                        #print(lf[0][:-1],where_at[:-1],lf[1][:-1],line[:-1],count_l,count_k)
                        #print('0000'*50)
                        #exit()
                        print(leaf1[:-1],leaf2[:-1] ,count_l,count_k,('-------------------------------'))
                        new_count_sibling +=  1.5* (count_l+count_k)
                        new_count_test +=  1.5* (count_l+count_k)
                        new_count_uncle +=  1.5* (count_l+count_k)
                        

                    taxa_mrca =current_node.taxa
                    lineage_middle_mrca = [ined[:-1] for ined in lineage_middle_mrca]
                    taxa_mrca =taxa_mrca.union(set(lineage_middle_mrca))
                    addr_mrca =essential.current_address(taxa_mrca,sp)
                    leaves_from_above =find_level_from_above(addr_mrca)
                    
                    expansion_factor = max(0,leaves_from_above-2)
                    lineage_middle =lineage_middle_full[:-1]
                    
                    mid_len= len(lineage_middle_mrca)
                    mid_len=0

    


                    flag=True
                    index=0
                    if flag_4:
                        if NNI_SP_test>=1 or hop+1>2:
                            lineage_middle =[]
                        else:
                            lineage_middle = lineage_middle_full[mid_len:-1]
                    
                    print('===>',len(lineage_middle),index)
                    print(lineage_middle)
                    while flag and index<len(lineage_middle):
                        line =lineage_middle[index]
                        index=index+1
                        print(leaf1[:-1],leaf2[:-1],where_at[:-1],line, ('-------------------------------'))
                        if line[-1]==';':
                            count_l,count_k= correct_accounting_by_product(data,leaf1[:-1],where_at[:-1],leaf2[:-1],line[:-1],sp_string)
                            #print(lf[0][:-1],where_at[:-1],lf[1][:-1],line[:-1],count_l,count_k)
                            #print('0000'*50)
                            #exit()
                        else:
                            count_l,count_k= correct_accounting_by_product(data,leaf1[:-1],where_at[:-1],leaf2[:-1],line,sp_string)

                        
                        if (count_l+count_k) <2:
                            flag=False

                        else:
                            if not essential.isequal_set(line,test_comparision_sibling) and not essential.isequal_set(line,test_comparison_uncle):
                                new_count_test +=  1.5* (count_l+count_k)
                                
                            if sibling_test_count:
                                if not essential.isequal_set(line,sibling_comparision_sibling):
                                    new_count_sibling +=  1.5* (count_l+count_k)
                                
                            if uncle_test_count:        
                                if not essential.isequal_set(line,uncle_comparision_sibling):
                                    new_count_uncle +=  1.5* (count_l+count_k)
                            
                    '''
                    for line in lineage_middle:
                        #sibling_list= essential.find_sibling(sp,[])
                        #uncle_list= essential.find_uncle(sp,[])
                        if not essential.isequal_set(line,comparison_sibling) and not essential.isequal_set(line,comparison_uncle):
                            count_l,count_k= correct_accounting_by_product(data,leaf1[:-1],where_at[:-1],leaf2[:-1],line[:-1],sp_string)
                            #print(lf[0][:-1],where_at[:-1],lf[1][:-1],line[:-1],count_l,count_k)
                            #print('0000'*50)
                            #exit()
                            print(leaf1[:-1],leaf2[:-1] ,count_l,count_k,('-------------------------------'))
                            
                            if (count_l+count_k) <10:
                                break
                            else:
                                new_count +=  1.5* (count_l+count_k)
                        
                        if hop==0:
                            new_couasasant +=  scale* (count_l+count_k)
                        else:
                            new_count +=  scale* (count_l+count_k)

                    '''
    print('367'*10,new_count_test,new_count_sibling,new_count_uncle)        
    return new_count_test,new_count_sibling,new_count_uncle


'''
def recurssive_correction(old_clade,where_at,scale):
    new_count=0
    for lf  in leaves:
        tr1= red.parse(lf)
        tr1.label_internal()
        if tr1.isLeaf:
            continue
        else:
            lfe= extract_leaves(lf)
        if  not essential.isequal(lf,old_clade) and not essential.isequal(lf,leaves[0]) and not essential.isequal(lf,leaves[1]):
            for line in lineage_middle:
                count_l,count_k= correct_accounting_by_product(data,lfe[0][:-1],where_at[:-1],lfe[1][:-1],line[:-1],sp_string)
                #print(lf[0][:-1],where_at[:-1],lf[1][:-1],line[:-1],count_l,count_k)
                #print('0000'*50)
                #exit()
                new_count +=  scale* (count_l+count_k)
            #exit()
            
'''        

def search_address(sorted_grouped, old_clade, where_at,flag):
    
    df_filtered_where_at = sorted_grouped[
        sorted_grouped[flag].apply(lambda x: essential.isequal(x, where_at))
    ]
    pprint.pprint(df_filtered_where_at.to_dict())
    to_return = {'uncle': [], 'sib': [], 'uncle_idx': [], 'sib_idx': []}
    
    if flag=='What_moved':
        old_clade=essential.map_pair_sibling(where_at,sibling_list)
        
    for idx, row in df_filtered_where_at.iterrows():
        com_sib = row['comparison_sibling']
        com_uncle = row['comparison_uncle']

        if pd.notna(com_sib) and essential.isequal(old_clade, com_sib):
            to_return['sib'].append(row)
            to_return['sib_idx'].append(idx)

        if pd.notna(com_uncle) and essential.isequal(old_clade, com_uncle):
            to_return['uncle'].append(row)
            to_return['uncle_idx'].append(idx)

    print(to_return)
    return to_return
    pass

def search_address_reverse(sorted_grouped, old_clade, where_at,flag):
    
    df_filtered_where_at = sorted_grouped[
        sorted_grouped['Where_at'].apply(lambda x: essential.isequal(x, where_at))
    ]
    pprint.pprint(df_filtered_where_at.to_dict())
    to_return = {'uncle': [], 'sib': [], 'uncle_idx': [], 'sib_idx': []}
    
    old_clade_1=essential.map_pair_sibling(old_clade,sibling_list)
        
    for idx, row in df_filtered_where_at.iterrows():
        com_sib = row['comparison_sibling']
        com_uncle = row['comparison_uncle']
        What_moved = row['What_moved']

        if pd.notna(com_sib) and essential.isequal(old_clade, com_sib):
            to_return['sib'].append(row)
            to_return['sib_idx'].append(idx)

        if pd.notna(com_uncle) and essential.isequal(old_clade, com_uncle):
            to_return['uncle'].append(row)
            to_return['uncle_idx'].append(idx)

    print(to_return)
    return to_return
    pass



def correct_count(data,sorted_grouped,sp_string):
    
    flag = 'What_moved'
    flag_f='Where_at'
    
    #print(sorted_grouped[flag])
    sorted_grouped['Sib_leaf_'+flag] = sorted_grouped[flag].apply(lambda x: is_leaf(x))
    sorted_grouped['Sib_leaf_'+flag_f] = sorted_grouped[flag_f].apply(lambda x: is_leaf(x))
    print(sorted_grouped['Sib_leaf_'+flag])
    
    #sorted_grouped= correct_count(data,sorted_grouped,'What_moved','total_count',sp_string,)
    #sorted_grouped= correct_count(data,sorted_grouped,'comparison_sibling','sibling_count',sp_string)
    #sorted_grouped= correct_count(data,sorted_grouped,'comparison_uncle','uncle_count',sp_string)
    
    #sorted_grouped['Sib_leaf'] = sorted_grouped['What_moved'].apply(lambda x: is_leaf(x))
    #leaf_grouped = (sorted_grouped.groupby('Where_at'))
    #print(sorted_grouped)
    sorted_grouped['total_count'] = sorted_grouped['total_count'].astype(float)
    sorted_grouped['sibling_count'] = sorted_grouped['sibling_count'].astype(float)
    sorted_grouped['uncle_count'] = sorted_grouped['uncle_count'].astype(float)
    sp= red.parse(sp_string)
    sp.label_internal()
    
    for idx, row in sorted_grouped.iterrows():
        if not row['Sib_leaf_'+flag] and row[flag]:
            old_count= row['total_count']
            old_clade =row[flag]
            where_at = row['Where_at']
            NNI_SP_test =row['NNI_sp']
            NNI_SP_sib,NNI_SP_uncle =None, None
            if NNI_SP_test==0:
                continue
            row_querry = search_address(sorted_grouped, old_clade, where_at,'Where_at')
            #NNI_SP_test,NNI_SP_sib,NNI_SP_uncle
            row_uncle = row_sib = None
            idx_uncle = idx_sib = None
            sibling_test_count = uncle_test_count = None

            if len(row_querry['uncle']) > 0:
                row_uncle = row_querry['uncle'][0]
                idx_uncle = row_querry['uncle_idx']

                uncle_test_count = row_uncle['What_moved']
                old_count_uncle = row_uncle['uncle_count']
                NNI_SP_uncle =row_uncle['NNI_sp']

            if len(row_querry['sib']) > 0:
                row_sib = row_querry['sib'][0]
                idx_sib = row_querry['sib_idx']
                

                sibling_test_count = row_sib['What_moved']
                old_count_sib = row_sib['sibling_count']
                NNI_SP_sib =row_sib['NNI_sp']

            print('rok', row_uncle, row_sib)
            print('rok', idx_uncle, idx_sib)

            #exit()
            
            test_count = [row['comparison_sibling'],row['comparison_uncle']]
            
            
            
            print(old_clade)
            old_1 =red.parse(old_clade)
            old_1.label_internal()
            curr=  essential.current_address(old_1.taxa,sp)
            if curr.parent:
                sib_1= [aer  for aer in curr.parent.children if aer !=curr][0]
                uncle=sib_1.to_newick()
            leaves= extract_leaves(old_clade)
            new_count =0
            #for leaf  in leaves:
            print('---'*40,where_at,old_clade,NNI_SP_test,leaves)
            
            count_x,count_y,count_z,taxa_mrca= correct_accounting(data,leaves[0][:-1],where_at[:-1],leaves[1][:-1],uncle[:-1])
            print(count_x,count_y,count_z)
            
            lineage_middle =extract_mrca_linage(sp_string,old_clade,where_at[:-1])
            lineage_middle = [ined[:-1] for ined in lineage_middle]
            taxa_mrca =taxa_mrca.union(set(lineage_middle))
            addr_mrca =curr.parent#essential.current_address(taxa_mrca,sp)
            #print(taxa_mrca)
            #print(addr_mrca)
            scale=1#.5+ (0.5)*(find_level_from_above(addr_mrca)+find_level_from_below(addr_mrca)+len(lineage_middle)-1)
            #print(scale,find_level_from_above(addr_mrca),find_level_from_below(addr_mrca))
            #exit()
            #if NNI_SP==1:
            #scale=2.20727*NNI_SP*(len(old_1.taxa)-1)*(len(lineage_middle)-1)
            ##print('e',count_x,count_y,count_z,find_level_from_above(addr_mrca),find_level_from_below(addr_mrca),len(lineage_middle),NNI_SP,'target')
            new_count += scale* (count_x+count_y)
            print('*****'*40,old_clade,where_at[:-1],lineage_middle)
            #sorted_grouped.loc[idx,count_flag]=old_count+new_count
            #lfe =  extract_all_lineages(old_clade)
            
            #print(lfe,old_clade)
            
            new_count_uncle,new_count_sibling,new_count_test=0,0,0
            if NNI_SP_test>=1:
                new_count_test,new_count_sibling,new_count_uncle =recurssive_correction(data,old_1,where_at,scale,lineage_middle,sp_string,NNI_SP_test,NNI_SP_sib,NNI_SP_uncle,test_count,sibling_test_count,uncle_test_count,False)
            
            sorted_grouped.loc[idx,'total_count']=old_count+1.5*new_count+new_count_test
                
                
            if len(row_querry['uncle'])>0:
                u_uncle=old_count_uncle+1.5*new_count+new_count_uncle
                for ied in idx_uncle:
                    sorted_grouped.loc[ied,'uncle_count']=u_uncle
            
            if len(row_querry['sib'])>0:
                u_sib= old_count_sib+1.5*new_count+new_count_sibling
                for ied in idx_sib:
                    sorted_grouped.loc[ied,'sibling_count']=u_sib
                
    print('84356'*100)
    for idx, row in sorted_grouped.iterrows():
        if not row['Sib_leaf_'+flag_f] and row[flag_f]:
            old_count= row['total_count']
            old_clade =row[flag_f]
            where_at = row[flag]
            NNI_SP_test =row['NNI_sp']
            NNI_SP_sib,NNI_SP_uncle =None, None
            if NNI_SP_test==0:
                continue
            #row_querry = search_address_reverse(sorted_grouped, old_clade, where_at,flag)
            row_querry = search_address_reverse(sorted_grouped, where_at,old_clade,'Where_at')
            #NNI_SP_test,NNI_SP_sib,NNI_SP_uncle
            print(row_querry,flag,old_clade, where_at)
            row_uncle = row_sib = None
            idx_uncle = idx_sib = None
            sibling_test_count = uncle_test_count = None

            if len(row_querry['uncle']) > 0:
                row_uncle = row_querry['uncle'][0]
                idx_uncle = row_querry['uncle_idx']

                uncle_test_count = row_uncle['What_moved']
                old_count_uncle = row_uncle['uncle_count']
                NNI_SP_uncle =row_uncle['NNI_sp']

            if len(row_querry['sib']) > 0:
                row_sib = row_querry['sib'][0]
                idx_sib = row_querry['sib_idx']
                

                sibling_test_count = row_sib['What_moved']
                old_count_sib = row_sib['sibling_count']
                NNI_SP_sib =row_sib['NNI_sp']

            print('rok', row_uncle, row_sib)
            print('rok', idx_uncle, idx_sib)

            #exit()
            
            test_count = [row['comparison_sibling'],row['comparison_uncle']]
            print(test_count)
            print('12234'*855)
            print(row)
            print(idx)
            
            print(old_clade)
            old_1 =red.parse(old_clade)
            old_1.label_internal()
            curr=  essential.current_address(old_1.taxa,sp)
            if curr.parent:
                sib_1= [aer  for aer in curr.parent.children if aer !=curr][0]
                uncle=sib_1.to_newick()
            leaves= extract_leaves(old_clade)
            new_count =0
            #for leaf  in leaves:
            print('---'*40,where_at,old_clade,NNI_SP_test,leaves)
            
            count_x,count_y,count_z,taxa_mrca= correct_accounting(data,leaves[0][:-1],where_at[:-1],leaves[1][:-1],uncle[:-1])
            print(count_x,count_y,count_z)
            
            lineage_middle =extract_mrca_linage(sp_string,old_clade,where_at[:-1])
            lineage_middle = [ined[:-1] for ined in lineage_middle]
            taxa_mrca =taxa_mrca.union(set(lineage_middle))
            addr_mrca =curr.parent#essential.current_address(taxa_mrca,sp)
            #print(taxa_mrca)
            #print(addr_mrca)
            scale=1#.5+ (0.5)*(find_level_from_above(addr_mrca)+find_level_from_below(addr_mrca)+len(lineage_middle)-1)
            #print(scale,find_level_from_above(addr_mrca),find_level_from_below(addr_mrca))
            #exit()
            #if NNI_SP==1:
            #scale=2.20727*NNI_SP*(len(old_1.taxa)-1)*(len(lineage_middle)-1)
            ##print('e',count_x,count_y,count_z,find_level_from_above(addr_mrca),find_level_from_below(addr_mrca),len(lineage_middle),NNI_SP,'target')
            new_count += scale* (count_x+count_y)
            print('*****'*40,old_clade,where_at[:-1],lineage_middle)
            #sorted_grouped.loc[idx,count_flag]=old_count+new_count
            #lfe =  extract_all_lineages(old_clade)
            
            #print(lfe,old_clade)
            new_count_uncle,new_count_sibling,new_count_test=0,0,0
            
            #if NNI_SP_test>=1 and row['Sib_leaf_'+flag]:
            new_count_test,new_count_sibling,new_count_uncle =recurssive_correction(data,old_1,where_at,scale,lineage_middle,sp_string,NNI_SP_test,NNI_SP_sib,NNI_SP_uncle,test_count,sibling_test_count,uncle_test_count,True)
           
            
            sorted_grouped.loc[idx,'total_count']=old_count+1.5*new_count+new_count_test
                
                
            if len(row_querry['uncle'])>0:
                u_uncle=old_count_uncle+1.5*new_count+new_count_uncle
                for ied in idx_uncle:
                    sorted_grouped.loc[ied,'uncle_count']=u_uncle
            
            if len(row_querry['sib'])>0:
                u_sib= old_count_sib+1.5*new_count+new_count_sibling
                for ied in idx_sib:
                    sorted_grouped.loc[ied,'sibling_count']=u_sib
                

    return sorted_grouped
                    
            
            

         
    




def parse1():
    parser = argparse.ArgumentParser(description="IQTree on Simphy and dupcoal")
    parser.add_argument('--sp', type=str, help="Species tree")
    parser.add_argument('--gt', type=str, help="Gene tree list")
    #parser.add_argument('--lineages',type=lambda s: s.split('/'),help="'/'-separated list of lineages (e.g. l1/l2)")
    parser.add_argument('--output', type=str, help="Name of output file")
    parser.add_argument('--direction', type=int,default=0, help="Run DAFT_Direction.py")
    parser.add_argument('--excel', type=int,default=0, help="Produce Excel")
    parser.add_argument('--sibling', type=str, default='0', help="Run sibling test")
    
    args = parser.parse_args()
    return args


parser = parse1()
essential= daft_essential()
reco =reconcils()
red= readWrite.readWrite()
Il=ILS.ILS()

sp_string = parser.sp
#lineages = parser.lineages 
gene_treefile =parser.gt
produce_excel=parser.excel
run_direction=parser.direction
#print(id_it(sp_string))
#exit()
#lineage1= lineages[0]
#lineage2= lineages[1]
sibling_flag=parser.sibling
out_file=parser.output#+'_'+lineage1+'_'+lineage2


lis=[]
siblings = []
what_moved = []
total_count_=[]


data=pd.read_csv(gene_treefile, sep=',').to_numpy()

sp= red.parse(sp_string)
sp.label_internal()

list_sp_tree=reco.get_get_current_lineages(sp,[])
newick_lineage= address_to_newick(list_sp_tree)

sib_lineage= find_sib_lineage_pair(data)


result = accounting(data,sib_lineage)    
sibling_list= essential.find_sibling(sp,[])
uncle_list= essential.find_uncle(sp,[])


grouped,sorted_grouped,df = sorting_arrangement()


sorted_grouped= put_sibling_uncle(grouped,sorted_grouped,df)
sorted_grouped=put_z(sorted_grouped,grouped)
node_map,branch_map,labeled_sp= essential.id_it(sp_string)
sorted_grouped_converted= essential.idfy_it(sorted_grouped,node_map)


pd.DataFrame(branch_map).to_csv('branch_map.csv',index=False)
#print(sorted_grouped)
#exit()
write_significance(sorted_grouped_converted,out_file,labeled_sp)

if produce_excel:
    convert_excel(out_file)
if run_direction:
    call_direction(sorted_grouped,gene_treefile,sp,out_file)
        
clean_folder(out_file)

    
