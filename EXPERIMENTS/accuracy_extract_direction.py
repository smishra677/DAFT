import re
from pathlib import Path

import pandas as pd


def daft_folder_relations_to_newick(
    folder,
    branch_map_file="branch_map.csv",
    output_file="Summary_direction.csv",
    start_idx=1,
    end_idx=10,
    file_prefix="DAFT_Direction_",
):
    """
    Read DAFT_Direction_1.txt ... DAFT_Direction_10.txt from a folder.

    For each existing file:
        - extract RECEIVER / Donor pairs
        - if BIDIRECTIONAL, add both directions
        - add idx column from the file number
        - convert IDs back to Newick using branch_map.csv

    Missing files get one empty row with only idx present.

    Returns
    -------
    pandas.DataFrame with columns:
        Receiver_ID
        Donor_ID
        Receiver_Lineage
        Donor_Lineage
        idx
    """

    folder = Path(folder)

    # -----------------------------
    # Load branch_map.csv
    # -----------------------------
    branch_df = pd.read_csv(folder / branch_map_file)
    branch_df = branch_df[["id", "To"]].copy()

    # Important: use pandas nullable integer type
    # This allows merging with columns that may contain pd.NA
    branch_df["id"] = branch_df["id"].astype("Int64")

    # -----------------------------
    # Regex patterns
    # -----------------------------
    block_pattern = re.compile(
        r"INFERRED\s+RELATIONS\s*\(ONLY\s+NNI\s*>\s*1\)\s*:\s*=+\s*(.*?)\*+",
        flags=re.S | re.I
    )

    relation_pattern = re.compile(
        r"^\s*\d+\s+AND\s+\d+\s+RECEIVER\s*:\s*(\d+)\s+AND\s+Donor\s*:\s*(\d+)(.*)$",
        flags=re.I
    )

    rows = []
    # -----------------------------
    # Loop over DAFT_Direction_1.txt ... DAFT_Direction_10.txt
    # -----------------------------
    for idx in range(start_idx, end_idx + 1):
        daft_file = folder / f"{file_prefix}{idx}.txt"
        #print(daft_file,daft_file.exists())
        #exit()
        if not daft_file.exists():
            rows.append(
                {
                    "idx": idx,
                    "Receiver_ID": pd.NA,
                    "Donor_ID": pd.NA,
                }
            )

            continue

        with open(daft_file, "r") as f:
            daft_text = f.read()
            #print(daft_text)
            #exit()

        block_match = block_pattern.search(daft_text)

        if not block_match:
            rows.append(
                {
                    "idx": idx,
                    "Receiver_ID": pd.NA,
                    "Donor_ID": pd.NA,
                }
            )
            continue

        block = block_match.group(1)

        found_relation = False

        for line in block.splitlines():
            match = relation_pattern.search(line)

            if not match:
                continue

            found_relation = True

            receiver_id = int(match.group(1))
            donor_id = int(match.group(2))
            extra_text = match.group(3)

            # receiver first, donor second
            rows.append(
                {
                    "idx": idx,
                    "Receiver_ID": receiver_id,
                    "Donor_ID": donor_id,
                }
            )

            # if bidirectional, also add donor first, receiver second
            if "(BIDIRECTIONAL)" in extra_text.upper():
                rows.append(
                    {
                        "idx": idx,
                        "Receiver_ID": donor_id,
                        "Donor_ID": receiver_id,
                    }
                )

        if not found_relation:
            rows.append(
                {
                    "idx": idx,
                    "Receiver_ID": pd.NA,
                    "Donor_ID": pd.NA,
                }
            )

    # -----------------------------
    # Build relations DataFrame
    # -----------------------------
    relations_df = pd.DataFrame(rows)

    if relations_df.empty:
        result_df = pd.DataFrame(
            columns=[
                "Receiver_ID",
                "Donor_ID",
                "Receiver_Lineage",
                "Donor_Lineage",
                "idx",
            ]
        )

        if output_file is not None:
            result_df.to_csv(folder / output_file, index=False)

        return result_df

    # Important: make these nullable Int64, not object
    relations_df["Receiver_ID"] = relations_df["Receiver_ID"].astype("Int64")
    relations_df["Donor_ID"] = relations_df["Donor_ID"].astype("Int64")

    # -----------------------------
    # Add Receiver_Lineage
    # -----------------------------
    result_df = relations_df.merge(
        branch_df.rename(
            columns={
                "id": "Receiver_ID",
                "To": "Receiver_Lineage",
            }
        ),
        on="Receiver_ID",
        how="left"
    )

    # -----------------------------
    # Add Donor_Lineage
    # -----------------------------
    result_df = result_df.merge(
        branch_df.rename(
            columns={
                "id": "Donor_ID",
                "To": "Donor_Lineage",
            }
        ),
        on="Donor_ID",
        how="left"
    )

    # -----------------------------
    # Check missing IDs in branch_map.csv
    # Ignore intentionally empty rows from missing files
    # -----------------------------
    missing_receiver = result_df[
        result_df["Receiver_ID"].notna()
        & result_df["Receiver_Lineage"].isna()
    ]

    missing_donor = result_df[
        result_df["Donor_ID"].notna()
        & result_df["Donor_Lineage"].isna()
    ]

    if not missing_receiver.empty:
        missing_ids = missing_receiver["Receiver_ID"].unique().tolist()
        raise KeyError(f"Receiver IDs not found in branch_map.csv: {missing_ids}")

    if not missing_donor.empty:
        missing_ids = missing_donor["Donor_ID"].unique().tolist()
        raise KeyError(f"Donor IDs not found in branch_map.csv: {missing_ids}")

    # -----------------------------
    # Final column order, idx last
    # -----------------------------
    result_df = result_df[
        [
            "Receiver_ID",
            "Donor_ID",
            "Receiver_Lineage",
            "Donor_Lineage",
            "idx",
        ]
    ]

    result_df = result_df.sort_values(
        by=["idx", "Receiver_ID", "Donor_ID"],
        na_position="last"
    ).reset_index(drop=True)

    if output_file is not None:
        result_df.to_csv(folder / output_file, index=False)

    return result_df

A=None

B=None

C=None


val_key = [
(A,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_M',[('A;','F;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_N',[('A;','H;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_O',[('H;','A;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_P',[('C;','H;')]),
(A,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_Q',[('A;','F;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_R',[('A;','H;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_S',[('H;','A;')]),
(A,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_Q_4',[('A;','F;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_R_4',[('A;','H;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_S_4',[('H;','A;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_T',[('C;','H;')]),
(A,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_Q_6',[('A;','F;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_R_6',[('A;','H;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_S_6',[('H;','A;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_T_5',[('C;','H;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_T_6',[('C;','H;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_T_7',[('C;','R;')]),
(A,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_U',[('A;','F;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_V',[('A;','H;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_W',[('H;','A;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_X',[('C;','H;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_Y',[('C;','F;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_Z',[('C;','M;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_AA',[('A;','H;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_BB',[('A;','H;')]),
(A,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_CC',[('A;','F;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_DD',[('A;','H;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_EE',[('H;','A;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_FF',[('C;','H;')]),
(A,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_CC_1',[('A;','F;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_DD_1',[('A;','H;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_EE_1',[('H;','A;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_FF_1',[('C;','H;')]),
(A,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_GG',[('C;','(G,F);')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_HH',[('H;','(A,B);')]),
(A,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_II',[('(G,F);','C;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_JJ',[('(A,B);','H;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_KK',[('(M,N);','(I,C);')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_LL',[('(O,((Q,R),P));','((E,(D,H)),((C,I),(G,F)));')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_MM',[('P;','T;'),('E;','(C,I);')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_NN',[('C;','H;'),('Q;','S;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_OO',[('J;','G;'),('K;','G;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_PP',[('Q;','O;'),('Q;','T;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_QQ',[('C;','H;'),('((E,(D,H)),((C,I),(G,F)));','(L,K);')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_QQ_1',[('(V,W);','(N,M);'),('Q;','W;')]),
(A,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_RR',[('A;','F;')]),
(B,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_SS',[('A;','H;')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_TT',[('(M,N);','(I,C);')]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_Gante_corrected',[]),
(C,'/N/slate/samishr/daft-e/RESULTS/RESULTS/MAY_19/DAFT_results_Gante_uncorrected',[])
]


for key,val,add in val_key:
    print(val)
    try:
        daft_folder_relations_to_newick(val)
    except:
        continue