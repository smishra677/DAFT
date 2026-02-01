import numpy as np
import os, sys
sys.path.append("./DAFT_utils/reconcILS")
sys.path.append("./DAFT_utils/")
from reconcILS import *
from utils_reconcILS import *
import pandas as pd 
import argparse
from DAFT_essential import *

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
        #print(gt)
        gt =gt[0]
        #print(gt)
        gt =gt.replace('e-', '0')
        gt_tr= red.parse(gt)
        gt_tr.label_internal()


        if  essential.current_address(taxa_1,gt_tr) and essential.current_address(taxa_2,gt_tr):
            filtered_data.append([gt])
    
    return filtered_data


def call_direction(sorted_grouped,gene_treefile,sp):
    list_unique= filter_direction(sorted_grouped,-2.25,10,0)
    #list_unique=pairs
    sibling_list= essential.find_sibling(sp,[])

    list_non_unique=add_sibling_bidirection(sorted_grouped,list_unique,sibling_list)

    if len(list_unique)==0:
        print('Zero Significant')

    else:
        #list_non_unique=pairs1

        script = "./DAFT_Direction.py"
        sp = sp_string
        output = "out"

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

def accounting(data):
    result={'Where_at':[],'What_moved':[],'NNI_sp':[],'total_count':[]}
    sib_lineage= find_sib_lineage_pair(data)
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
            #result1['NNI_sp']+=[dist]
            #result1['total_count']+=[valu]

            result1['Where_at']+=[key1_t.to_newick()]
            result1['What_moved']+=[key_t.to_newick()]
            #result1['NNI_sp']+=[dist]
            #result1['total_count']+=[valu]
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


def put_sibling_uncle(grouped,sorted_grouped,df):
    sorted_grouped['comparison_sibling'] = sorted_grouped['What_moved'].apply(lambda x: essential.map_pair(x,sibling_list))
    sorted_grouped['comparison_uncle'] = sorted_grouped['What_moved'].apply(lambda x: essential.map_pair(x,uncle_list))


    for where_at, group in sorted_grouped.groupby('Where_at',as_index=False):
        for idx, row in group.iterrows():
            What_moved=row['What_moved']
            comp_uncle= row['comparison_uncle']
            comp_sibling=row['comparison_sibling']
            NNI_sp =row['NNI_sp']
            #print(where_at,comp_uncle,What_moved)
            if NNI_sp==0:
                sorted_grouped.loc[idx, 'comparison_sibling'] = None
                sorted_grouped.loc[idx, 'comparison_uncle'] = None
                continue
            
            if essential.isequal_set(comp_sibling,where_at) or essential.is_in(where_at,comp_sibling):
                sorted_grouped.loc[idx, 'comparison_sibling'] = None

            if essential.isequal_set(comp_uncle,where_at) or essential.is_in(where_at,comp_uncle):
                sorted_grouped.loc[idx, 'comparison_uncle'] = None
                continue
                
            #print(What_moved,comp_uncle)
            if comp_uncle:
                Tree1= red.parse(What_moved)
                Tree2= red.parse(comp_uncle)
                Tree3= red.parse(where_at)

                Tree1.label_internal()
                Tree2.label_internal()
                Tree3.label_internal()

                distance1=abs(essential.findDist(sp,Tree1.taxa,Tree3.taxa)-2)
                distance2=abs(essential.findDist(sp,Tree2.taxa,Tree3.taxa)-2)

                #print('None',What_moved,comp_uncle,where_at,distance1,distance2)
                if distance2>distance1:
                    sorted_grouped.loc[idx, 'comparison_uncle'] = None

    return sorted_grouped

def filter_direction(sorted_grouped,Z,tc,nni):
    filtered_df_1 = sorted_grouped[
        ((sorted_grouped['Z-value-uncle'] <= Z) | (sorted_grouped['Z-value-sibling'] <= Z))&(sorted_grouped['total_count']>tc)&(sorted_grouped['NNI_sp']>nni)
    ]

    pairs = list(zip(filtered_df_1['What_moved'], filtered_df_1['Where_at']))
    return pairs


def put_z(sorted_grouped,grouped):
    sorted_grouped['sibling_count'] = (
        sorted_grouped.apply(essential._sibling_count, axis=1,args=(grouped,)).astype('Int64')
    )

    essential.compute_z(sorted_grouped,'sibling',df)
    sorted_grouped['uncle_count'] = (
        sorted_grouped.apply(essential._uncle_count, axis=1,args=(grouped,)).astype('Int64')
    )



    essential.compute_z(sorted_grouped,'uncle',df)
    sorted_grouped = (
        sorted_grouped
        .groupby('Where_at', group_keys=False, dropna=False)
        .apply(essential._fix_group_nested)
    )

    return sorted_grouped


def write_significance(sorted_grouped):
    sorted_grouped = sorted_grouped.sort_values(by=[ 'NNI_sp'])

    cols = [
        'NNI_sp',
        'Test_lineage',
        'total_count',
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
            'total_count',
            'comparison_uncle',
            'uncle_count',
            'Z-value-uncle'
        ]

    sorted_grouped.rename(columns={"What_moved": "Test_lineage"}, inplace=True)


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
    SPECIAL_AFTER = ['total_count', 'Z-value-uncle']  # put 3 tabs after these columns
    TAB_STR = '\t' * 5
    SPACE_SEP = ' ' * padding
    TABSTOP = 8  

    with open('./DAFT_Significance.txt', 'w') as oe:
        oe.write(f"Species Tree = {sp_string}\n")
        oe.write("=" * 40 + "\n")

        for where_at, group in sorted_grouped.groupby('Where_at', as_index=False):
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


def convert_excel():
        script = "./DAFT_produce_excel.py"
        argv = [
                sys.executable,         
                script,
                
            ]
        status = os.spawnv(os.P_WAIT, sys.executable, argv)
        exit_code = os.WEXITSTATUS(status) if hasattr(os, "WEXITSTATUS") else (status >> 8)
        if exit_code != 0:
            raise RuntimeError(f"DAFT_produce_excel.py failed with code {exit_code}")

def clean_folder():
        script = "./clean_folder.py"
        argv = [
                sys.executable,         
                script,
                
            ]
        status = os.spawnv(os.P_WAIT, sys.executable, argv)
        exit_code = os.WEXITSTATUS(status) if hasattr(os, "WEXITSTATUS") else (status >> 8)
        if exit_code != 0:
            raise RuntimeError(f"clean_folder.py failed with code {exit_code}")
            
          
    




def parse1():
    parser = argparse.ArgumentParser(description="IQTree on Simphy and dupcoal")
    parser.add_argument('--sp', type=str, help="Species tree")
    parser.add_argument('--gt', type=str, help="Gene tree list")
    parser.add_argument('--excel', type=int,default=0, help="Produce Excel")
    parser.add_argument(
        '--lineages',
        type=lambda s: s.split('/'),
        help="Comma-separated list of lineages (e.g. l1,l2)"
    )
    parser.add_argument('--output', type=str, help="Name of output file")
    parser.add_argument('--sibling', type=str, default='0', help="Species tree")
    
    args = parser.parse_args()
    return args


parser = parse1()
essential= daft_essential()
reco =reconcils()
red= readWrite.readWrite()
Il=ILS.ILS()

sp_string = parser.sp
lineages = parser.lineages 
gene_treefile =parser.gt
produce_excel=parser.excel

lineage1= lineages[0]
lineage2= lineages[1]
sibling_flag=parser.sibling
out_file=parser.output+'_'+lineage1+'_'+lineage2

lis=[]
siblings = []
what_moved = []
total_count_=[]
#newick_lineage=[]

data=pd.read_csv(gene_treefile, sep=',').to_numpy()

sp= red.parse(sp_string)
sp.label_internal()



data=filter_data(data,lineage1,lineage2)
list_sp_tree=reco.get_get_current_lineages(sp,[])
newick_lineage= address_to_newick(list_sp_tree)


result = accounting(data)    
sibling_list= essential.find_sibling(sp,[])
uncle_list= essential.find_uncle(sp,[])


grouped,sorted_grouped,df = sorting_arrangement()
sorted_grouped= put_sibling_uncle(grouped,sorted_grouped,df)
sorted_grouped=put_z(sorted_grouped,grouped)

write_significance(sorted_grouped)

if produce_excel:
    convert_excel()
call_direction(sorted_grouped,gene_treefile,sp)
        
clean_folder()

    