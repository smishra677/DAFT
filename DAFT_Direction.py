import os
import sys
import pandas as pd
import ast
import argparse
import numpy as np
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

def read_list_file(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        return ast.literal_eval(f.read())
    
    
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


def parse1():
    parser = argparse.ArgumentParser(description="Input to DAFT Direction")

    parser.add_argument('--sp', type=str, default=None, help="Species tree as Newick text")
    parser.add_argument('--sp_file', type=str, default=None, help="File containing species tree")
    parser.add_argument('--gt', type=str, help="List of gene trees")
    parser.add_argument('--path', type=str, default="./DAFT_utils", help="Path to DAFT_utils")
    parser.add_argument('--verbose', type=int, default=1, help="Verbose mode. 1 = yes, 0 = no")

    parser.add_argument(
        '--lineages_file',
        type=str,
        required=True,
        help="File containing list of lineage tuples"
    )

    parser.add_argument(
        '--lineagesN_file',
        type=str,
        default=None,
        help="Optional file containing non-unique lineage tuples"
    )

    parser.add_argument('--output', type=str, help="Name of output file")
    
    parser.add_argument('--rooting', type=str, default=None, help="Outgroup taxon or comma-separated outgroup taxa")
    parser.add_argument('--forced', type=int, default=0, help="Force rerooting using --rooting")
    parser.add_argument('--allow_inconsistent_rooting', type=int, default=0, help="Continue even if rooting is inconsistent")
    parser.add_argument('--cache_hash', type=str, default=None, help="Input hash used for djiNNI cache")
    parser.add_argument('--ignore_duplication', type=int, default=0, help="Ignore gene tree containing duplication and continue")
    parser.add_argument('--random_seed', type=int, default=42, help="Random seed for djiNNI")
    
    
    args = parser.parse_args()

    if args.sp_file:
        with open(args.sp_file, "r", encoding="utf-8") as f:
            args.sp = f.read().strip()

    if not args.sp:
        parser.error("Please provide species tree using --sp or --sp_file")

    if not args.sp.endswith(";"):
        args.sp += ";"

    args.lineages = read_list_file(args.lineages_file)

    if args.lineagesN_file:
        args.lineagesN = read_list_file(args.lineagesN_file)
    else:
        args.lineagesN = []

    return args

parser = parse1()

path = parser.path

sys.path.append(path)
sys.path.append(path + "/reconcILS")

from reconcILS import *
from utils_reconcILS import *
from DAFT_essential import *
from DAFT_validate import *




total_count_cutoff = 5
indiv_count_cutoff = 5

reco =reconcils()
red= readWrite.readWrite()
essential= daft_essential()
validate= daft_validate()

sp_string = parser.sp
gene_treefile = essential.convert_gt_to_csv(parser.gt)
lineages = parser.lineages
lineages_bidrectional = parser.lineagesN
out_filec = parser.output
verbose=parser.verbose

ignore_duplicate =parser.ignore_duplication
force_root=parser.forced
allow_incos=parser.allow_inconsistent_rooting
to_root=parser.rooting
random_seed=parser.random_seed
np.random.seed(int(random_seed))


        
if parser.cache_hash:
    djiNNI_hash = parser.cache_hash
    validation = validate.validateme('djiNNI',sp_string,gene_treefile,out_filec,parser,to_root,force_root,allow_incos,ignore_duplicate,verbose,hashed=1)

    
else:
    validation = validate.validateme('djiNNI',sp_string,gene_treefile,out_filec,parser,to_root,force_root,allow_incos,ignore_duplicate,verbose,hashed=0)

    sp_string = validation.species_tree_newick
    gene_treefile = validation.gene_treefile
    djiNNI_hash = validation.input_hash





def convert_species(df,sp_string):
    '''
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
    '''
    species_to_letters  = {
    "D_albomicans": "A",
    "D_nasuta": "B",
    "D_kepulauana": "C",
    "D_neonasuta": "D",
    "D_sulfurigaster_albostrigata": "E",
    "D_pulaua": "F",
    "D_sulfurigaster_bilimbata": "G",
    "D_sulfurigaster_sulfurigaster": "H",
    "D_neohypocausta": "I",
    "D_immigrans_kari17": "K",
    "D_immigrans": "J",
    "D_pruinosa": "L",
    "D_arawakana": "M",
    "D_dunni": "N",
    "D_cardini": "O",
    "D_ornatifrons": "P",
    "D_subbadia": "Q",
    "D_pallidipennis": "R",
    "D_funebris": "S",
    "D_guttifera": "T",
    "D_innubila": "U",
    "D_mush_saotome": "V",
    "D_quadrilineata": "W"  }
    
    #letter_to_species ={value:key for key,value in species_to_letters.keys()}
    sp_= red.parse(sp_string)
    sp_string= sp_.to_newick_change(species_to_letters)
    for idx, row in df.iterrows():

        # ---- Primary fields ----
        old_lineage1 = row['Lineage1']
        old_lineage2 = row['Lineage2']
        old_what_moved = row['What_moved']
        old_to_where = row['To_where']
        
        print(old_lineage1,old_lineage2)

        # ---- Minor sibling fields ----
        old_minor_moved = row['Minor_Moved']
        old_minor_sibling = row['minor_sibling']
        old_minor_sib_what = row['minor_sibling_What_moved']
        old_minor_sib_to_where = row['minor_sibling_To_where']

        # ---- Significant pairs ----
        old_sig = row['Significant_Pairs']
        old_receiver, old_donor = '', ''

        if pd.notna(old_sig) and 'AND' in old_sig:
            parts = old_sig.split('AND')
            old_receiver = parts[0].strip()
            old_donor = parts[1].strip()

        # ---- Helper ----
        def convert_tree(val):
            if pd.notna(val) and val != '':
                tree = red.parse(val)
                return tree.to_newick_change(species_to_letters)
            return val

        # ---- Convert primary ----
        new_lineage1 = convert_tree(old_lineage1)
        new_lineage2 = convert_tree(old_lineage2)
        new_what_moved = convert_tree(old_what_moved)
        new_to_where = convert_tree(old_to_where)

        # ---- Convert significant ----
        new_receiver = convert_tree(old_receiver)
        new_donor = convert_tree(old_donor)

        # ---- Convert minor ----
        new_minor_moved = convert_tree(old_minor_moved)
        new_minor_sibling = convert_tree(old_minor_sibling)
        new_minor_sib_what = convert_tree(old_minor_sib_what)
        new_minor_sib_to_where = convert_tree(old_minor_sib_to_where)

        # ---- Store back (ALWAYS write back) ----
        df.at[idx, 'Lineage1'] = new_lineage1
        df.at[idx, 'Lineage2'] = new_lineage2
        df.at[idx, 'What_moved'] = new_what_moved
        df.at[idx, 'To_where'] = new_to_where

        df.at[idx, 'Minor_Moved'] = new_minor_moved
        df.at[idx, 'minor_sibling'] = new_minor_sibling
        df.at[idx, 'minor_sibling_What_moved'] = new_minor_sib_what
        df.at[idx, 'minor_sibling_To_where'] = new_minor_sib_to_where

        df.at[idx, 'Significant_Pairs'] = str(new_receiver) + ' AND ' + str(new_donor)

    return df, sp_string

def write_direction(df,network_output,sp,out_filec):
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
    with open("DAFT_Direction_"+out_filec+".txt", "w") as f:
        f.write("## DAFT Direction\n\n")
        
        f.write("SIGNIFICANT PAIRS\n")
        f.write("=" * 80 + "\n")
        for i1, i2,NNI_ in df[['Lineage1', 'Lineage2','NNI_']].itertuples(index=False, name=None):
            #f.write(f"BETWEEN {i1} AND {i2}  {' '*(25-len(i1)+len(i2))+str(essential.find_dist_string(sp,i1,i2))} NNI APART \n")
            f.write(f"BETWEEN {i1} AND {i2}  {' '*(25-len(i1)+len(i2))+str(NNI_)} NNI APART \n")
        f.write("*" * 80 + "\n")
        f.write("\n") 
        f.write("\n") 
        f.write("\n")
        

        # Write Data Table 1
        f.write("DATA TABLE1 (ONLY NNI > 1) \n")
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
            NNI_between = row['NNI_']
            if NNI_between>1:
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
        f.write("DATA TABLE2 (BIDIRECTIONAL) (ONLY NNI > 1) \n")
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
            NNI_between = row['NNI_']
            if NNI_between>1:
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
        f.write("INFERRED RELATIONS (ONLY NNI > 1) :\n")
        f.write("=" * 80 + "\n")
        for ide, row in df.iterrows():
            NNI_between = row['NNI_']
            if NNI_between>1:
                try:
                    z = row.get("Z_score", np.nan)
                    zs = row.get("Z_score_sibling", np.nan)
                    zf = float(z)
                    zsf = float(zs)
                    is_bidirectional = (not pd.isna(zf)) and (not pd.isna(zsf)) and (zsf<-1.96)
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
        f.write('Network (ONLY NNI > 1) :  ' + network_output+'\n')
        f.write("*" * 80 + "\n")
        f.write("\n") 
        f.write("\n") 
        f.write("\n")


def aggregrate_lineage(dfe,lineage1,lineage2,output):
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
    aggregated_data.to_csv(f'introgression_1_group_{output}.csv', index=True)

    df = pd.read_csv(f"./introgression_1_group_{output}.csv", sep=',')
    return df,count_lineage1,count_lineage2

def collapse_clade(df):
    groups = []
    used = set()
    #print(df)

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
    #print(df)
    #print(lineage1,lineage2)
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
           
def find_equal(visited,querry):
    for tup in visited:
        if tuple_equal(tup, querry):
            return True
    return False
        
        

def run_tranform(data,out_filec,lineages_bidrectional,sp_string,list_df,path,verbose,djiNNI_hash,random_seed):
    vis=[]
    for lineageM, lineage2 in lineages_bidrectional:
        #print('Testing Direction for: ',lineageM,lineage2)
        for lineage1 in lineageM:
            if find_equal(vis,(lineage1,lineage2)):
                continue
            vis.append((lineage1,lineage2))
            l1 = red.parse(lineage1)
            l2= red.parse(lineage2)
            
            l1.label_internal()
            l2.label_internal()
            taxa_1=l1.taxa
            taxa_2=l2.taxa
            filtered_data=[]
            index_data=[]

            for idex, gt in enumerate(data):
                gt =gt[0]
                gt_tr= red.parse(gt)
                gt_tr.label_internal()
                if  essential.current_address(taxa_1,gt_tr) and essential.current_address(taxa_2,gt_tr) and reco.get_current_sister(gt_tr,taxa_1)==taxa_2:
                    filtered_data.append(gt)
                    index_data.append(idex)
                
            
                    


            if len(filtered_data)!=0:
                

                script = "./DAFT_Transform.py"  
                #sp = sp_string
                lineages1 = f"{lineage1}/{lineage2}"   
                #print(lineages1)

                #print('asasas',verbose)
                #output = out_filec
                djiNNI_output_name= essential.dji_output_name(out_filec,sp_string,lineage1,lineage2)
                
                argv = [
                    sys.executable, script,
                    "--sp", str(sp_string),
                    "--lineages", str(lineages1),
                    "--path", path,
                    "--verbose", str(verbose),
                    "--gt_list", *filtered_data,
                    "--gt_list_index", *map(str, index_data),
                    "--cache_hash", str(djiNNI_hash),
                    "--random_seed", str(random_seed),
                    "--output", djiNNI_output_name,
                ]
                status = os.spawnv(os.P_WAIT, sys.executable, argv)
                exit_code = status if os.name == "nt" else (status >> 8) 

                if exit_code != 0:
                    raise RuntimeError(f"DAFT_Transform.py exited with code {exit_code}")

                df=pd.read_csv('./djiNNI_'+djiNNI_output_name+'.csv',sep=',')
                #print(df)
                
                #import pprint
                #pprint.pprint(df.to_dict()['Sibling'])
                #print('-->')
                #pprint.pprint(df.to_dict()['What_moved'])
                #print(sp_tree_lineages)
                
                df["Sibling"] = df["Sibling"].apply(lambda v: essential.match_lineage(v, sp_tree_lineages))
                df["What_moved"] = df["What_moved"].apply(lambda v: essential.match_lineage(v, sp_tree_lineages))
                
                #total_count_=[]
                #newick_lineage=[]
                #print(df)
                #print(df)

                df=df.dropna()
                dfe =df.to_numpy()
                sp= red.parse(sp_string)
                sp.label_internal()
                
                
                df,count_lineage1,count_lineage2= aggregrate_lineage(dfe,lineage1,lineage2,out_filec)
                #print(list_df)
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



#print(lineages)
#exit()

#print(lineages_bidrectional)

#data=pd.read_csv(gene_treefile, sep=',').to_numpy()
data = validate.read_gene_tree_data(gene_treefile)


sp =red.parse(sp_string)
sp.label_internal()

list_df={'Labeled_species':[],'Total_gene_trees':[],'Lineage1':[],'Count1':[],'Lineage2':[],'Count2':[],'What_moved':[],'To_where':[],'Minor_Moved':[],'Minor_moved_count':[]}
sp_tree_lineages= essential.find_all_lineage(sp)

if len(lineages_bidrectional)==0:
    lineages_bidrectional= add_sibling_bidirection(lineages,sp)


#print(lineages_bidrectional)
run_tranform(data,out_filec,lineages_bidrectional,sp_string,list_df,path,verbose,djiNNI_hash,random_seed)

df = pd.DataFrame(list_df)
#df.to_csv(f'results1_{out_filec}.csv',index=False)



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

#df,sp_string1= convert_species(df,sp_string)
#df.to_csv(f'results1_{out_filec}.csv',index=False)

df['NNI_'] = df.apply(
    lambda row: essential.find_dist_string(sp, row['Lineage1'], row['Lineage2']),
    axis=1
)


node_map,branch_map,sp_labeled= essential.id_it(sp_string)

pd.DataFrame(branch_map).to_csv('branch_map.csv',index=False)
df_converted= essential.idfy_it_direction(df,node_map)
#print(sp_labeled.to_newick())
#exit()
#put it in a network
df_net = df[(df['NNI_']>1)]

sp_labeled =essential.put_network_in_tree(df_net,sp_labeled)
network_output=essential.to_network(sp_labeled)

#print(df)fp
#print(sp_labeled.to_newick())
write_direction(df_converted,network_output,sp,out_filec)
clean_folder(out_filec)