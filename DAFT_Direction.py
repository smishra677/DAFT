import os
import sys
sys.path.append("./DAFT_utils")
sys.path.append("./DAFT_utils/reconcILS")
from reconcILS import *
from utils_reconcILS import *
import pandas as pd 
import ast
from DAFT_essential import *
import argparse
import numpy as np



reco =reconcils()
red= readWrite.readWrite()
essential= daft_essential()



def write_direction(df,network_output,sp):
    padding = 2
    SPACE_SEP = " " * padding
    TAB_STR = "\t" * 5
    TABSTOP = 8
    
    SPECIAL_AFTER = ["To_where", "minor_sibling_To_where"]

    TABLE1_COL =[
        'Significant_Pairs',
        'Total_gene_trees',
        'Lineage1',
        'Count1',
        'Lineage2',
        'Count2',
        'Major_Moved'
    ]

    TABLE2_COL =[
    'Significant_Pairs',
    'Minor_Moved',
    'Minor_sibling',
    'Total_gene_trees',
    'Minor_sibling_count',
    'Minor_moved_count',
    'Z_score'
    ]

    conversion_dic1 ={
        'Significant_Pairs':'Significant_Pairs',
        'Total_gene_trees':'Total_gene_trees',
        'Lineage1':'Lineage1',
        'Count1':'Count1',
        'Lineage2':'Lineage2',
        'Count2':'Count2',
        'Major_Moved':'What_moved'

    }


    conversion_dic2 ={
        'Significant_Pairs':'Significant_Pairs',
        'Minor_Moved':'Minor_Moved',
        'Minor_sibling':'minor_sibling',
        'Total_gene_trees':'minor_sibling_Total_gene_trees',
        'Minor_sibling_count':'minor_sibling_count',
        'Minor_moved_count':'Minor_moved_count',
        'Z_score':'Z_score_sibling'

    }

    widths_base = essential.compute_widths(df, TABLE1_COL)
    widths_full = essential.compute_widths(df, TABLE2_COL)


    df = df.fillna("-")
    df.replace("", "-", inplace=True)
    #print(df)
    #df.to_csv('test.csv',index=False)
    #exit()
    with open("DAFT_Direction.txt", "w") as f:
        f.write("## DAFT Direction\n\n")
        
        f.write("SIGNIFICANT PAIRS\n")
        f.write("=" * 80 + "\n")
        for i1, i2 in df[['Lineage1', 'Lineage2']].itertuples(index=False, name=None):
            f.write(f"BETWEEN {i1} AND {i2}  {' '*(25-len(i1)+len(i2))+str(essential.find_dist_string(sp,i1,i2))} NNI APART \n")
        f.write("*" * 80 + "\n")
        f.write("\n") 
        f.write("\n") 
        f.write("\n")
        

        # Write Data Table 1
        f.write("DATA TABLE1\n")
        #print([TABLE1_COL,widths_full])
        
        
        header = ""
        
        for i, c in enumerate(TABLE1_COL, 1):
            header += f"{c:<{widths_base[c]}}"
            if i < len(TABLE1_COL):
                header += TAB_STR if c in SPECIAL_AFTER else SPACE_SEP
        f.write("=" * len(header.expandtabs(TABSTOP).rstrip()) + "\n")
        f.write(header+'\n')
        f.write("=" * len(header.expandtabs(TABSTOP).rstrip()) + "\n")

        #subset=df 
        
        for ide, row in df.iterrows():
            line = ""
            for i, c in enumerate(TABLE1_COL, 1):
                val = row[conversion_dic1[c]] if conversion_dic1[c] in row else np.nan
                s = essential.render_val(c, val)
                line += f"{s:<{widths_base[c]}}"
                if i < len(TABLE1_COL):
                    line += TAB_STR if c in SPECIAL_AFTER else SPACE_SEP
            line += "\n"
            f.write(line)
        f.write("*" * len(header.expandtabs(TABSTOP).rstrip()) + "\n")
        
        f.write("\n") 
        f.write("\n") 
        f.write("\n")  

        
        # Write Data Table 2 (Bidirectional Data)
        f.write("DATA TABLE2 (BIDIRECTIONAL)\n")
        header = ""
        for i, c in enumerate(TABLE2_COL, 1):
            header += f"{c:<{widths_full[c]}}"
            if i < len(TABLE2_COL):
                header += TAB_STR if c in SPECIAL_AFTER else SPACE_SEP
        f.write("=" * len(header.expandtabs(TABSTOP).rstrip()) + "\n")
        f.write(header+'\n')
        f.write("=" * len(header.expandtabs(TABSTOP).rstrip()) + "\n")

        #subset=df 
        for ide, row in df.iterrows():
            line = ""
            for i, c in enumerate(TABLE2_COL, 1):
                val = row[conversion_dic2[c]] if conversion_dic2[c] in row else np.nan
                s = essential.render_val(c, val)
                line += f"{s:<{widths_full[c]}}"
                if i < len(TABLE2_COL):
                    line += TAB_STR if c in SPECIAL_AFTER else SPACE_SEP
            line += "\n"
            f.write(line)
        
        f.write("*" * len(header.expandtabs(TABSTOP).rstrip()) + "\n")
        f.write("\n") 
        f.write("\n") 
        f.write("\n")

        # Write Data Table 2 (Bidirectional Data)
        f.write("INFERRED RELATIONS:\n")
        f.write("=" * 80 + "\n")
        for ide, row in df.iterrows():
            try:
                z = row.get("Z_score", np.nan)
                zs = row.get("Z_score_sibling", np.nan)
                zf = float(z)
                zsf = float(zs)
                is_bidirectional = (not pd.isna(zf)) and (not pd.isna(zsf)) and (zf > zsf)
            except Exception:
                is_bidirectional = False

            if is_bidirectional:
                line = f"{row['Significant_Pairs']:<25}{' RECEIVER:'+row['What_moved']+' AND Donor:'+ row['Minor_Moved'] +'  (BIDIRECTIONAL)'}\n"
            else:
                line = f"{row['Significant_Pairs']:<25}{' RECEIVER:'+row['What_moved']+' AND Donor:'+ row['Minor_Moved']}\n"
            f.write(line)

        f.write("*" * 80 + "\n")
        f.write("\n") 
        f.write("\n")
        f.write("=" * 80 + "\n")
        f.write('Network: ' + network_output+'\n')
        f.write("*" * 80 + "\n")
        f.write("\n") 
        f.write("\n") 
        f.write("\n")


def aggregrate_lineage(dfe,lineage1,lineage2):
    count_lineage1=0
    count_lineage2=0
    leni=[]

    for entry in dfe:
        #print(entry[-2],lineage1)
        if (essential.isequal(lineage1,entry[-2]) and essential.isequal(lineage2,entry[-3])):
            leni.append(entry)
            count_lineage1=count_lineage1+1

        if (essential.isequal(lineage2,entry[-2]) and essential.isequal(lineage1,entry[-3])):
            leni.append(entry)
            count_lineage2=count_lineage2+1


    df= pd.DataFrame(leni,columns=['idx','Replicate','Path','From_Where_moved','Sibling','What_moved','NNI'])

    aggregated_data = df.groupby(['From_Where_moved','What_moved','Sibling']).agg(
        total_count=('Replicate', 'count')
    )
    aggregated_data.to_csv('introgression_1_group.csv', index=True)

    df = pd.read_csv("./introgression_1_group.csv", sep=',')
    return df,count_lineage1,count_lineage2

def collapse_clade(df):
    groups = []
    used = set()

    for i, row_i in df.iterrows():
        if i in used:
            continue
        current_group = [i]
        for j, row_j in df.iterrows():
            if j <= i or j in used:
                continue
            if (essential.isequal_set(row_i["What_moved"], row_j["What_moved"]) and
                essential.isequal_set(row_i["Sibling"], row_j["Sibling"])):
                current_group.append(j)
        groups.append(current_group)
        used.update(current_group)

    # Collapse groups
    collapsed = []
    for g in groups:
        sub = df.loc[g]
        collapsed.append({
            "From_Where_moved": sub["From_Where_moved"].iloc[0], 
            "What_moved": sub["What_moved"].iloc[0],
            "Sibling": sub["Sibling"].iloc[0],
            "total_count": sub["total_count"].sum()
        })

    df = pd.DataFrame(collapsed)
    return df

def get_donor_receipient(lineage1,lineage2,count_lineage1,count_lineage2,choice):
    rename_map = {
        lineage1[:-1]: 'receiver',
        lineage2[:-1]: 'donor'
    }
    if count_lineage1>count_lineage2:
        rename_map = {
        lineage1[:-1]: 'receiver',
        lineage2[:-1]: 'donor'}
    elif count_lineage1<count_lineage2:
        rename_map = {
        lineage1[:-1]: 'donor',
        lineage2[:-1]: 'receiver'}
    else:
        if choice==0:
            rename_map = {
            lineage1[:-1]: 'receiver',
            lineage2[:-1]: 'donor'}
        else:
            rename_map = {
            lineage1[:-1]: 'donor',
            lineage2[:-1]: 'receiver'}
                
            
    return rename_map



def update_big_list(df,data,list_df,lineage1,lineage2,rename_map,choice):
    count_nni_score ={1:0,2:0,3:'',4:''}
    lineage_seen=[]
    for idx,line in enumerate(df.to_numpy()):
        count_nni_score[idx+3]=line[1]
        count_nni_score[idx+1]=line[3]
        lineage_seen.append(line[1])

    list_df['To_where']+=[line[0]]
    list_df['Total_gene_trees']+=[len(data)]

    if count_nni_score[3]=='':
        if lineage1 not in lineage_seen:
            count_nni_score[3]=lineage1
        else:
            count_nni_score[3]=lineage2

    if count_nni_score[4]=='':
        if lineage1 not in lineage_seen:
            count_nni_score[4]=lineage1
        else:
            count_nni_score[4]=lineage2

    list_df['Lineage1']+=[count_nni_score[3]]
    list_df['Lineage2']+=[count_nni_score[4]]
    list_df['Count1']+=[count_nni_score[1]]
    list_df['Count2']+=[count_nni_score[2]]




    updated_newick = essential.rename_subtrees(sp_string, rename_map)


    list_df['Labeled_species']+=[updated_newick]
    if count_nni_score[1]>count_nni_score[2]:
        list_df['What_moved']+=[list_df['Lineage1'][-1]]
        list_df['Minor_Moved']+=[list_df['Lineage2'][-1]]
        list_df['Minor_moved_count']+=[count_nni_score[2]]
    elif count_nni_score[1]<count_nni_score[2]:
        list_df['What_moved']+=[list_df['Lineage2'][-1]]
        list_df['Minor_Moved']+=[list_df['Lineage1'][-1]]
        list_df['Minor_moved_count']+=[count_nni_score[1]]
    else:
        if choice==0:
            list_df['What_moved']+=[list_df['Lineage1'][-1]]
            list_df['Minor_Moved']+=[list_df['Lineage2'][-1]]
            list_df['Minor_moved_count']+=[count_nni_score[2]]
        else:
            list_df['What_moved']+=[list_df['Lineage2'][-1]]
            list_df['Minor_Moved']+=[list_df['Lineage1'][-1]]
            list_df['Minor_moved_count']+=[count_nni_score[1]]
            
            

def run_tranform(lineages_bidrectional,sp_string,list_df):
    for lineageM, lineage2 in lineages_bidrectional:
        #print('Testing Direction for: ',lineageM,lineage2)
        for lineage1 in lineageM:
            l1 = red.parse(lineage1)
            l2= red.parse(lineage2)
            l1.label_internal()
            l2.label_internal()
            taxa_1=l1.taxa
            taxa_2=l2.taxa
            filtered_data=[]

            for gt in data:
                gt =gt[0]
                gt_tr= red.parse(gt)
                gt_tr.label_internal()
                if  essential.current_address(taxa_1,gt_tr) and essential.current_address(taxa_2,gt_tr) and reco.get_current_sister(gt_tr,taxa_1)==taxa_2:
                    filtered_data.append(gt)


            if len(filtered_data)!=0:

                script = "./DAFT_Transform.py"  
                #sp = sp_string
                lineages1 = f"{lineage1}/{lineage2}"   



                output = "out"
                argv = [
                    sys.executable, script,
                    "--sp", str(sp_string),
                    "--lineages", str(lineages1),
                    "--gt_list", *filtered_data,
                    "--output", str('out'),
                ]
                status = os.spawnv(os.P_WAIT, sys.executable, argv)
                exit_code = status if os.name == "nt" else (status >> 8) 
                if exit_code != 0:
                    raise RuntimeError(f"DAFT_Transform.py exited with code {exit_code}")


                df=pd.read_csv('./introgression_out.csv',sep=',')
                
                df["Sibling"] = df["Sibling"].apply(lambda v: essential.match_lineage(v, sp_tree_lineages))
                df["What_moved"] = df["What_moved"].apply(lambda v: essential.match_lineage(v, sp_tree_lineages))
                
                #total_count_=[]
                #newick_lineage=[]
                

                df=df.dropna()
                dfe =df.to_numpy()
                sp= red.parse(sp_string)
                sp.label_internal()
                
                df,count_lineage1,count_lineage2= aggregrate_lineage(dfe,lineage1,lineage2)
                #print(df)
                df= collapse_clade(df)


                # ressolve equal count for donor and receipient  
                choice= np.random.choice([0,1])
                rename_map= get_donor_receipient(lineage1,lineage2,count_lineage1,count_lineage2,choice)
                
                update_big_list(df,filtered_data,list_df,lineage1,lineage2,rename_map,choice)


                

            else:
                list_df['Labeled_species']+=[sp_string]
                list_df['Total_gene_trees']+=[0]
                list_df['Lineage1']+=[lineage1]
                list_df['Lineage2']+=[lineage2]
                list_df['Count1']+=[0]
                list_df['Count2']+=[0]
                list_df['What_moved']+=[None]
                list_df['To_where']+=[None]
                list_df['Minor_Moved']+=[None]
                list_df['Minor_moved_count']+=[0]
                
                #list_df={'Labeled_species':[],'Total_gene_trees':[],'Lineage1':[],'Count1':[],'Lineage2':[],'Count2':[],'What_moved':[],'To_where':[]}



            #print(list_df)
            #exit()
    #return df


def add_sibling_bidirection(pairs,sp):
    sibling_list= essential.find_sibling(sp,[])
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



def parse1():
    parser = argparse.ArgumentParser(description="Input to DAFT Direction")
    parser.add_argument('--sp', type=str, help="Species tree")
    parser.add_argument('--gt', type=str, help="List of gene trees")
    parser.add_argument(
        '--lineages',
        type=lambda s: ast.literal_eval(s),
        help="List of lineage tuples to check direction, e.g. \"[(1,2), (3,4)]\""
    )
    parser.add_argument(
        '--lineagesN', default='[]',
        type=lambda s: ast.literal_eval(s),
        help="List of lineage tuples, e.g. \"[(1,2), (3,4)]\""
    )
    parser.add_argument('--output', type=str, help="Name of output file")
    
    args = parser.parse_args()
    return args


parser = parse1()
sp_string = parser.sp
gene_treefile =parser.gt
lineages = parser.lineages
lineages_bidrectional = parser.lineagesN
out_filec='big_output'
total_count_cutoff=15
indiv_count_cutoff=0

#print(lineages)
#exit()

data=pd.read_csv(gene_treefile, sep=',').to_numpy()



sp =red.parse(sp_string)
list_df={'Labeled_species':[],'Total_gene_trees':[],'Lineage1':[],'Count1':[],'Lineage2':[],'Count2':[],'What_moved':[],'To_where':[],'Minor_Moved':[],'Minor_moved_count':[]}
sp_tree_lineages= essential.find_all_lineage(sp)

if len(lineages_bidrectional)==0:
    lineages_bidrectional= add_sibling_bidirection(lineages,sp)

run_tranform(lineages_bidrectional,sp_string,list_df)


df = pd.DataFrame(list_df)
df.to_csv('results1.csv',index=False)



#bidirectional processing
df = essential.sibling_bidirectional(df,lineages,sp_string)
df1 = pd.DataFrame(
    df[['Lineage1', 'Lineage2']].apply(lambda row: sorted([row['Lineage1'], row['Lineage2']]), axis=1).tolist(),
    columns=['Lineage1_s', 'Lineage2_s']
)
df1['original_index'] = df.index
unique_pairs = df1.drop_duplicates(subset=['Lineage1_s', 'Lineage2_s'])
df = df.loc[unique_pairs['original_index']]
df =df[df.inUnique==True]
df = df[(df['Count1']+df['Count2'])>total_count_cutoff &  (df['Count1']>indiv_count_cutoff) & (df['Count2']>indiv_count_cutoff)]


#put it in a network
sp_labeled =essential.put_network_in_tree(df,sp_string)
network_output=essential.to_network(sp_labeled)


write_direction(df,network_output,sp)