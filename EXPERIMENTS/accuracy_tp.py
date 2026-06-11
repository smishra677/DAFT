import sys
sys.path.append("../DAFT_utils/reconcILS/")
sys.path.append("../DAFT_utils/")
from DAFT_essential import *
import pandas as pd


essential= daft_essential()


def extract_summary(val, expected,flag):
    df = pd.read_csv(val + '/Summary.csv', header=0)
    if flag==1:
        df = df[df['flag'] == 1]
    else:
        df = df[df['flag'] == 2]

    big_update = []

    for a, b in expected:
        print(a, b)

        matched_idx = set()

        for _, row in df.iterrows():
            cond1 = essential.isequal(row['Focal_lineage'], a) and essential.isequal(row['Test_lineage'], b)
            cond2 = essential.isequal(row['Focal_lineage'], b) and essential.isequal(row['Test_lineage'], a)

            if cond1 or cond2:
                matched_idx.add(row['idx'])

        big_update.append(len(matched_idx))

    return sum(big_update), len(expected)*10


results_uncle={'Data_set':[],'Total_sims':[],'False_negatives':[],'True_positives':[],'False_negatives_rate':[],'True_positive_rate':[]}
results_sibling={'Data_set':[],'Total_sims':[],'False_negatives':[],'True_positives':[],'False_negatives_rate':[],'True_positive_rate':[]}



A='(((A,B),C),((F,G),H));'
B='(((((A,B),C),F),G),H);'
C=' (A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));'
#exit()
#print(extract_summary('./RESULTS/RESULTS/APRIL_19/DAFT_results_OO',[('H;','F;'),('(H,(G,F));','(A,B);')]))
#pprint.pprint(load_data_struct(list_sp_tree))
#exit()

val_key = [
(A,'./RESULTS/RESULTS/APRIL_19/DAFT_results_M',[('A;','F;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_N',[('A;','H;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_O',[('H;','A;')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_P',[('C;','H;')]),
(A,'./RESULTS/RESULTS/APRIL_19/DAFT_results_Q',[('A;','F;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_R',[('A;','H;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_S',[('H;','A;')]),
(A,'./RESULTS/RESULTS/APRIL_19/DAFT_results_Q_1',[('A;','F;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_R_1',[('A;','H;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_S_1',[('H;','A;')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_T',[('C;','H;')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_T_1',[('C;','H;')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_T_2',[('C;','R;')]),
(A,'./RESULTS/RESULTS/APRIL_19/DAFT_results_U',[('A;','F;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_V',[('A;','H;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_W',[('H;','A;')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_X',[('C;','H;')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_Y',[('C;','F;')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_Z',[('C;','M;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_AA',[('A;','H;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_BB',[('A;','H;')]),
(A,'./RESULTS/RESULTS/APRIL_19/DAFT_results_CC',[('A;','F;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_DD',[('A;','H;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_EE',[('H;','A;')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_FF',[('C;','H;')]),
(A,'./RESULTS/RESULTS/APRIL_19/DAFT_results_GG',[('C;','(G,F);')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_HH',[('H;','(A,B);')]),
(A,'./RESULTS/RESULTS/APRIL_19/DAFT_results_II',[('(G,F);','C;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_JJ',[('(A,B);','H;')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_KK',[('(M,N);','(I,C);')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_LL',[('(O,((Q,R),P));','((E,(D,H)),((C,I),(G,F)));')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_MM',[('P;','T;'),('E;','(C,I);')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_NN',[('C;','H;'),('Q;','S;')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_OO',[('J;','G;'),('K;','G;')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_PP',[('Q;','O;'),('Q;','T;')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_QQ',[('(V,W);','(N,M);'),('Q;','W;')]),
(A,'./RESULTS/RESULTS/APRIL_19/DAFT_results_RR',[('A;','F;')]),
(B,'./RESULTS/RESULTS/APRIL_19/DAFT_results_SS',[('A;','H;')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_TT',[('(M,N);','(I,C);')]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_Gante_corrected',[]),
(C,'./RESULTS/RESULTS/APRIL_19/DAFT_results_Gante_uncorrected',[])
]

for _sp_string,val,expected in val_key:
    tp,total_test= extract_summary(val,expected,1)
    #data_struct= extract_summary(val,expected)
    sim_id=val.split('DAFT_results_')[-1].strip()
    results_uncle['Data_set']+=[sim_id]
    results_uncle['Total_sims']+=[total_test]
    results_uncle['False_negatives']+=[total_test-tp]
    results_uncle['True_positives']+=[tp]
    results_uncle['False_negatives_rate']+=[(total_test-tp)/(total_test)]
    results_uncle['True_positive_rate']+=[tp/(total_test)]
    
    
    tp,total_test= extract_summary(val,expected,2)
    sim_id=val.split('DAFT_results_')[-1].strip()
    results_sibling['Data_set']+=[sim_id]
    results_sibling['Total_sims']+=[total_test]
    results_sibling['False_negatives']+=[total_test-tp]
    results_sibling['True_positives']+=[tp]
    results_sibling['False_negatives_rate']+=[(total_test-tp)/(total_test)]
    results_sibling['True_positive_rate']+=[tp/(total_test)]
    
pd.DataFrame(results_uncle).to_csv('Results_rate_uncle_M_TT.csv',index=False)
pd.DataFrame(results_sibling).to_csv('Results_rate_sibling_M_TT.csv',index=False)
