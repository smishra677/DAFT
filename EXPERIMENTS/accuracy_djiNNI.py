import sys
sys.path.append("./DAFT_utils/reconcILS/")
sys.path.append("./DAFT_utils/")
from DAFT_essential import *
import pandas as pd


essential= daft_essential()


def extract_summary_direction(val, receiver, donor):
    df = pd.read_csv(val + "/Summary_direction.csv", header=0)

    matched_idx = set()

    for _, row in df.iterrows():
        if pd.notna(row["Donor_Lineage"]) and pd.notna(row["Receiver_Lineage"]):

            cond1 = (
                essential.isequal(row["Donor_Lineage"], donor)
                and essential.isequal(row["Receiver_Lineage"], receiver)
            )

            if cond1:
                matched_idx.add(row["idx"])

    return len(matched_idx)


def extract_summary(results,val, expected,sim_id):
    # expected is a list of truth links in this order:
    #     (donor, receiver)

    df = pd.read_csv(val + '/Summary.csv', header=0)
    df = df[(df['flag'] == 1) | (df['flag'] == 2)]

    for donor, receiver in expected:
        matched_idx = set()

        for _, row in df.iterrows():
            focal = row['Focal_lineage']
            test = row['Test_lineage']

            cond1 = essential.isequal(focal, donor) and essential.isequal(test, receiver)
            cond2 = essential.isequal(focal, receiver) and essential.isequal(test, donor)

            if cond1 or cond2:
                matched_idx.add(row['idx'])

        detected = len(matched_idx)
        not_detected = 10 - detected

        if detected != 0:
            seen_direction = extract_summary_direction(val,receiver,donor)

            false_negative = detected - seen_direction
            true_positive = seen_direction
            false_negative_rate = false_negative / detected
            true_positive_rate = true_positive / detected

        else:
            false_negative = 0
            true_positive = 0
            false_negative_rate = 0
            true_positive_rate = 0

        results['Data_set'] += [sim_id]
        results['Receiver'] += [receiver]
        results['Donor'] += [donor]
        results['Detected_DAFT'] += [detected]
        results['Not_detected_DAFT'] += [not_detected]
        results['Total_test'] += [detected]
        results['False_negatives'] += [false_negative]
        results['True_positives'] += [true_positive]
        results['False_negatives_rate'] += [false_negative_rate]
        results['True_positive_rate'] += [true_positive_rate]

    return results


results={'Data_set':[],'Receiver':[],'Donor':[],'Detected_DAFT':[],'Not_detected_DAFT':[],'Total_test':[],'False_negatives':[],'True_positives':[],'False_negatives_rate':[],'True_positive_rate':[]}

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
    try:
        sim_id=val.split('DAFT_results_')[-1].strip()
        results= extract_summary(results,val,expected,sim_id)
    except Exception as err:
        print("Skipping", val, "because:", err)
        continue
    #data_struct= extract_summary(val,expected)
print(results)
pd.DataFrame(results).to_csv('Results_rate_djiNNI_M_TT.csv',index=False)
