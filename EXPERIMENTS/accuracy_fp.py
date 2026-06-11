import sys
sys.path.append("./DAFT_utils/reconcILS/")
sys.path.append("./DAFT_utils/")
from reconcILS import *
from utils_reconcILS import *
import pandas as pd
from DAFT_essential import *
import pprint
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm




essential= daft_essential()
reco =reconcils()
red= readWrite.readWrite()


def plotter(data_struct,val,flag):

    #D= ast.literal_eval(open('data_struct.txt').read())
    df = pd.DataFrame.from_dict(data_struct, orient="index")
    '''
    #df = pd.DataFrame(D)
    df = pd.DataFrame.from_dict(D, orient="index")
    # Replace empty strings with NaN
    #df = df.replace('', np.nan)

    # Ensure consistent row/column order
    df = df.reindex(index=df.columns, columns=df.columns)

    # Keep only upper triangular part
    df_upper = df.where(np.triu(np.ones(df.shape), k=0).astype(bool))

    print(df_upper)
    '''


    # =========================
    # 1. LOAD DATA
    # =========================

    #D = {...}  # your dictionary

    df_upper=df
    #df = df.iloc[1:, 1:]
    df.to_csv(val+'/summary_agg.csv')
    mask_x = df == ''
    df = df.replace('', 0).astype(float)
    Z = df.values

    valid_mask = ~mask_x.values
    Z_valid = Z[valid_mask]

    # =========================
    # 2. COMPUTE METRICS
    # =========================

    FP_count = np.sum(Z_valid[Z_valid > 0])
    #FN_count = -np.sum(Z_valid[Z_valid < 0])
    
    if flag==1:
        total_cells = np.sum(Z_valid >= 0)*10
    else:
        total_cells = np.sum(Z_valid >= 0)*5
    TP_cells = total_cells -FP_count

    #FN_cells = np.sum(Z_valid < 0)

    #total_cells = FP_count + FN_count + TP_cells
    #total_cells = FP_count + FN_count + TP_cells



    FP_rate = FP_count / total_cells if total_cells > 0 else 0
    #FN_rate = FN_count / total_cells if total_cells > 0 else 0
    TP_rate = TP_cells / total_cells if total_cells > 0 else 0

    # =========================
    # 3. PLOT
    # =========================
    mask_x = df_upper == ''
    df_upper = df_upper.replace('', 0).astype(float)
    Z = df_upper.values

    valid_mask = ~mask_x.values
    Z_valid = Z[valid_mask]
    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(1, 2, width_ratios=[3, 1])

    ax = fig.add_subplot(gs[0])

    max_abs = max(abs(Z_valid.min()), abs(Z_valid.max()))
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)

    x = np.arange(Z.shape[1] + 1)
    y = np.arange(Z.shape[0] + 1)

    m = ax.pcolormesh(
        x, y, Z,
        cmap="bwr",
        norm=norm,
        edgecolors="black",
        linewidth=0.6,
        shading="flat"
    )

    # Add X for impossible
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            if mask_x.iat[i, j]:
                ax.text(j + 0.5, i + 0.5, 'X',
                        ha='center', va='center',
                        color='black', fontsize=11, fontweight='bold')

    ax.set_xticks(np.arange(Z.shape[1]) + 0.5)
    ax.set_yticks(np.arange(Z.shape[0]) + 0.5)
    ax.set_xticklabels(df.columns, rotation=90)
    ax.set_yticklabels(df.index)
    ax.invert_yaxis()

    cbar = plt.colorbar(m, ax=ax)
    cbar.set_label("FP (>0 red) | TP (=0 white) | FN (<0 blue)")

    ax.set_title("Results(X- invalid introgression)")

    # =========================
    # 4. METRICS PANEL
    # =========================

    ax2 = fig.add_subplot(gs[1])
    ax2.axis("off")

    metrics_text = (
        f"Total Valid Cells acorss 10 runs : {total_cells}\n"
        f"True Positives (0) across 10 runs: {TP_cells}\n"
        f"Total FP Magnitude (>0) across 10 runs: {FP_count}\n\n"
        #f"Total FN Magnitude (<0) across 10 runs: {-FN_count}\n\n"
        f"TP Rate: {TP_rate:.6f}\n"
        f"FP Rate: {FP_rate:.6f}\n"
        #f"FN Rate: {-FN_rate:.3f}\n"
    )

    ax2.text(0.05, 0.95, metrics_text,
            va="top", fontsize=11)


    #plt.show()
    
    return total_cells,TP_cells,FP_count,TP_rate,FP_rate




def load_data_struct(list_sp_tree,flag):
    data_struct ={}
    for lin1 in list_sp_tree:
        if not  lin1.parent:
            continue
        lin1_newick=lin1.to_newick()
        lin1_sibling = [sib_ for sib_ in lin1.parent.children if sib_!=lin1][0]
        #list_sp_tree_newick+=[lin.to_newick()]
        data_struct[lin1_newick]={}
        for lin2 in list_sp_tree:
            if  not  lin2.parent:
                continue
            lin2_newick =lin2.to_newick()
            #if not essential.isequal(lin2.to_newick(), lin.to_newick()):
            if flag==2:
                if  not essential.is_in(lin1_sibling.to_newick(),lin2_newick) and not essential.is_in(lin2_newick,lin1_sibling.to_newick()) and not essential.is_in(lin1_newick, lin2_newick) and not essential.is_in(lin2_newick, lin1_newick) and not essential._is_sibling(lin1_newick,lin2_newick,sp_string) :
                    data_struct[lin1_newick][lin2_newick]=0
                    #data_struct[lin1.to_newick()][lin2_newick]=''
                    #data_struct[lin2_newick][lin1_newick]=''
                else:
                    data_struct[lin1_newick][lin2_newick]=''
                    #data_struct[lin1_newick][lin2_newick]=''
            else:
                if  not essential.is_in(lin1_sibling.to_newick(),lin2_newick) and not essential.is_in(lin2_newick,lin1_sibling.to_newick()) and not essential.is_in(lin1_newick, lin2_newick) and not essential.is_in(lin2_newick, lin1_newick) and not essential._is_uncle(lin1_newick,lin2_newick,sp_string) :
                    data_struct[lin1_newick][lin2_newick]=0
                    #data_struct[lin1.to_newick()][lin2_newick]=''
                    #data_struct[lin2_newick][lin1_newick]=''
                else:
                    data_struct[lin1_newick][lin2_newick]=''
                    #data_struct[lin1_newick][lin2_newick]=''
        
    return data_struct





def get_data(val,sp_string,flag):
    sp= red.parse(sp_string)

    sp.label_internal()
    df =pd.read_csv(val+'/Summary.csv',header=0)
    #print(df)
    
    if flag==1:
        df = df[df['flag'] == 1]
    else:
        df = df[df['flag'] == 2]
    list_sp_tree=reco.get_get_current_lineages(sp,[])
    #pprint.pprint(list_sp_tree)
    data_struct= load_data_struct(list_sp_tree,flag)
    #pprint.pprint(data_struct)

    #list_sp_tree_newick=[]

    #print(essential.is_in('A;','(A,B);'))
    #exit()
    super_data_struct=[]
    for _, group in df.groupby("idx"):
        #print(f"\nGroup idx = {value}")
        
        #if value in[1,3]:
        data= []
        for _, row in group.iterrows():
            focal = row["Focal_lineage"]
            test = row["Test_lineage"]
            #flag = row["flag"]
            
            data+=[(focal, test)]
        
        data_struct= load_data_struct(list_sp_tree,flag)
        for focal_com,test_com in data:
                #print(focal_com,test_com)
                if not pd.isna(focal_com):
                    
                    #print(data_struct)
                    key_focal = [ke for ke in  data_struct.keys() if essential.isequal(focal_com,ke)][0]
                    key_test = [ke for ke in data_struct.keys() if essential.isequal(test_com,ke)][0]
                    
                    #print(key_focal,key_test)
                    try:
                        data_struct[key_focal][key_test]+=1
                        continue
                    except:
                        Exception
                    try:
                        data_struct[key_test][key_focal]+=1
                        continue
                    except:
                        Exception
                    #print(data_struct[key_test][key_focal])
                    #exit()
                    
                    #super_data_struct+=[print(data_struct)]

                else:
                    continue
        super_data_struct+=[data_struct]
        
        #pprint.pprint(data_struct)

    data_struct= load_data_struct(list_sp_tree,flag)
    #print(data_struct)
    for struct in super_data_struct:
        for key1 in data_struct.keys():
            for key2 in data_struct.keys():
                #print(key1)
                data_struct[key1][key2]+= struct[key1][key2]
                #data_struct[key2][key1]+= struct[key2][key1]
                

    #pprint.pprint(data_struct)
    return data_struct
    
    
                



#pprint.pprint(load_data_struct(list_sp_tree))
#exit()

results_uncle={'Data_set':[],'Total_tests':[],'False_positives':[],'True_positives':[],'False_positive_rate':[],'True_positive_rate':[]}
results_sibling={'Data_set':[],'Total_tests':[],'False_positives':[],'True_positives':[],'False_positive_rate':[],'True_positive_rate':[]}



A='(((A,B),C),((F,G),H));'
B='(((((A,B),C),F),G),H);'
C=' (A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));'



val_key = [(A,'./RESULTS/RESULTS/APRIL_19/DAFT_results_A'),
           (B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_B'),
           (C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_C'),
           (A,'./RESULTS/RESULTS/APRIL_19/DAFT_results_D'),
           (B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_E'),
           (C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_F'),
           (A,'./RESULTS/RESULTS/APRIL_19/DAFT_results_G'),
           (B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_H'),
           (C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_I'),
           (A,'./RESULTS/RESULTS/APRIL_19/DAFT_results_J'),
           (B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_K'),
           (C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_L'),]



for sp_string,val in val_key:
    data_struct= get_data(val,sp_string,1)
    #exit()
    #print(data_struct)
    try:
        total_cells,TP_cells,FP_count,TP_rate,FP_rate = plotter(data_struct,val,1)
        print(total_cells,TP_cells,FP_count,TP_rate,FP_rate)
        sim_id=val.split('DAFT_results_')[-1].strip()
        results_uncle['Data_set']+=[sim_id]
        results_uncle['Total_tests']+=[total_cells]
        results_uncle['False_positives']+=[FP_count]
        results_uncle['True_positives']+=[TP_cells]
        results_uncle['False_positive_rate']+=[FP_rate]
        results_uncle['True_positive_rate']+=[TP_rate]
        #print(results)
    except:
        continue
        
#print(results)
pd.DataFrame(results_uncle).to_csv('Results_rate_uncle.csv',index=False)


for sp_string,val in val_key:
    data_struct= get_data(val,sp_string,2)
    #exit()
    #print(data_struct)
    try:
        total_cells,TP_cells,FP_count,TP_rate,FP_rate = plotter(data_struct,val,2)
        print(total_cells,TP_cells,FP_count,TP_rate,FP_rate)
        sim_id=val.split('DAFT_results_')[-1].strip()
        results_sibling['Data_set']+=[sim_id]
        results_sibling['Total_tests']+=[total_cells]
        results_sibling['False_positives']+=[FP_count]
        results_sibling['True_positives']+=[TP_cells]
        results_sibling['False_positive_rate']+=[FP_rate]
        results_sibling['True_positive_rate']+=[TP_rate]
    except:
        continue
        
pd.DataFrame(results_sibling).to_csv('Results_rate_sibling.csv',index=False)