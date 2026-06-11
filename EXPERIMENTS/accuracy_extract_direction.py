import re
from pathlib import Path

import pandas as pd


def parse_id_or_lineage(value):
    """
    If value is numeric, treat it as branch_map ID.
    Otherwise, treat it as already being the lineage/Newick string.
    """
    if pd.isna(value):
        return pd.NA, pd.NA

    value = str(value).strip().strip('"').strip("'").strip()

    if value == "":
        return pd.NA, pd.NA

    if re.fullmatch(r"\d+", value):
        return int(value), pd.NA

    return pd.NA, value


def daft_folder_relations_to_newick(folder,branch_map_file="DAFT_results_1/branch_map.csv", output_file="Summary_direction.csv", start_idx=1,end_idx=10,file_prefix="DAFT_Direction_",):
    """
    Read DAFT_Direction_1.txt ... DAFT_Direction_10.txt from a folder.

    Extract relations from lines like:

        A; AND F;    RECIPIENT:A; AND Donor:F;

    Also handles numeric branch IDs if they appear:

        1 AND 5      RECIPIENT:10 AND Donor:5;

    Rules
    -----
    - If RECIPIENT / Donor is numeric, convert using branch_map.csv.
    - If RECIPIENT / Donor is already lineage text like A;, F;, (G,F);,
      keep it directly.
    - If BIDIRECTIONAL appears, add both directions.
    """

    folder = Path(folder)

    # -----------------------------
    # Load branch_map.csv
    # -----------------------------
    branch_map_path = folder / branch_map_file

    if not branch_map_path.exists():
        raise FileNotFoundError(f"branch_map.csv not found: {branch_map_path}")

    branch_df = pd.read_csv(branch_map_path)
    branch_df = branch_df[["id", "To"]].copy()
    branch_df["id"] = pd.to_numeric(branch_df["id"], errors="coerce").astype("Int64")

    # -----------------------------
    # Regex patterns
    # -----------------------------
    block_pattern = re.compile(
        r"INFERRED\s+RELATIONS\s*\(ONLY\s+NNI\s*>\s*1\)\s*:\s*=+\s*(.*?)\*+",
        flags=re.S | re.I,
    )

    # This matches:
    # A; AND F;                 RECIPIENT:A; AND Donor:F;
    # 10 AND 5                  RECIPIENT:10 AND Donor:5
    #
    # It ignores the left side before RECIPIENT.
    # It extracts only the value after RECIPIENT: and after Donor:
    relation_pattern = re.compile(
        r"^\s*.*?\s+"
        r"(?:RECIPIENT|RECEPIENT)\s*:\s*(?P<receiver>.*?)\s+"
        r"AND\s+Donor\s*:\s*(?P<donor>.*?)\s*$",
        flags=re.I,
    )

    rows = []

    # -----------------------------
    # Loop over DAFT_Direction_1.txt ... DAFT_Direction_10.txt
    # -----------------------------
    for idx in range(start_idx, end_idx + 1):
        daft_file = folder / f"DAFT_results_{idx}/{file_prefix}{idx}.txt"

        if not daft_file.exists():
            rows.append(
                {
                    "idx": idx,
                    "Receiver_ID": pd.NA,
                    "Donor_ID": pd.NA,
                    "Receiver_Lineage": pd.NA,
                    "Donor_Lineage": pd.NA,
                }
            )
            continue

        with open(daft_file, "r", errors="ignore") as f:
            daft_text = f.read()

        block_match = block_pattern.search(daft_text)

        if not block_match:
            rows.append(
                {
                    "idx": idx,
                    "Receiver_ID": pd.NA,
                    "Donor_ID": pd.NA,
                    "Receiver_Lineage": pd.NA,
                    "Donor_Lineage": pd.NA,
                }
            )
            continue

        block = block_match.group(1)
        found_relation = False

        for line in block.splitlines():
            line = line.strip()

            if not line:
                continue

            # Example line:
            # A; AND F;                 RECIPIENT:A; AND Donor:F;
            #
            # Possible bidirectional line:
            # A; AND F;                 RECIPIENT:A; AND Donor:F; (BIDIRECTIONAL)

            is_bidirectional = bool(
                re.search(r"\bBIDIRECTIONAL\b", line, flags=re.I)
            )

            # Remove BIDIRECTIONAL text so donor is clean
            clean_line = re.sub(
                r"\s*\(?\s*BIDIRECTIONAL\s*\)?\s*$",
                "",
                line,
                flags=re.I,
            ).strip()

            match = relation_pattern.search(clean_line)

            if not match:
                continue

            found_relation = True

            receiver_raw = match.group("receiver").strip()
            donor_raw = match.group("donor").strip()

            receiver_id, receiver_lineage = parse_id_or_lineage(receiver_raw)
            donor_id, donor_lineage = parse_id_or_lineage(donor_raw)

            # Main direction
            rows.append(
                {
                    "idx": idx,
                    "Receiver_ID": receiver_id,
                    "Donor_ID": donor_id,
                    "Receiver_Lineage": receiver_lineage,
                    "Donor_Lineage": donor_lineage,
                }
            )

            # Reverse direction if BIDIRECTIONAL
            if is_bidirectional:
                rows.append(
                    {
                        "idx": idx,
                        "Receiver_ID": donor_id,
                        "Donor_ID": receiver_id,
                        "Receiver_Lineage": donor_lineage,
                        "Donor_Lineage": receiver_lineage,
                    }
                )

        if not found_relation:
            rows.append(
                {
                    "idx": idx,
                    "Receiver_ID": pd.NA,
                    "Donor_ID": pd.NA,
                    "Receiver_Lineage": pd.NA,
                    "Donor_Lineage": pd.NA,
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

    for col in [
        "Receiver_ID",
        "Donor_ID",
        "Receiver_Lineage",
        "Donor_Lineage",
        "idx",
    ]:
        if col not in relations_df.columns:
            relations_df[col] = pd.NA

    relations_df["Receiver_ID"] = pd.to_numeric(
        relations_df["Receiver_ID"], errors="coerce"
    ).astype("Int64")

    relations_df["Donor_ID"] = pd.to_numeric(
        relations_df["Donor_ID"], errors="coerce"
    ).astype("Int64")

    result_df = relations_df.copy()

    # -----------------------------
    # Fill Receiver_Lineage from branch_map only when Receiver_ID is numeric
    # -----------------------------
    receiver_map = branch_df.rename(
        columns={
            "id": "Receiver_ID",
            "To": "Receiver_Lineage_from_map",
        }
    )

    result_df = result_df.merge(
        receiver_map,
        on="Receiver_ID",
        how="left",
    )

    result_df["Receiver_Lineage"] = result_df["Receiver_Lineage"].fillna(
        result_df["Receiver_Lineage_from_map"]
    )

    result_df = result_df.drop(columns=["Receiver_Lineage_from_map"])

    # -----------------------------
    # Fill Donor_Lineage from branch_map only when Donor_ID is numeric
    # -----------------------------
    donor_map = branch_df.rename(
        columns={
            "id": "Donor_ID",
            "To": "Donor_Lineage_from_map",
        }
    )

    result_df = result_df.merge(
        donor_map,
        on="Donor_ID",
        how="left",
    )

    result_df["Donor_Lineage"] = result_df["Donor_Lineage"].fillna(
        result_df["Donor_Lineage_from_map"]
    )

    result_df = result_df.drop(columns=["Donor_Lineage_from_map"])

    # -----------------------------
    # Check missing numeric IDs only
    # Direct lineage strings do not need branch_map.csv
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
        missing_ids = missing_receiver["Receiver_ID"].dropna().unique().tolist()
        raise KeyError(f"Receiver IDs not found in branch_map.csv: {missing_ids}")

    if not missing_donor.empty:
        missing_ids = missing_donor["Donor_ID"].dropna().unique().tolist()
        raise KeyError(f"Donor IDs not found in branch_map.csv: {missing_ids}")

    # -----------------------------
    # Final column order
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
        by=[
            "idx",
            "Receiver_ID",
            "Donor_ID",
            "Receiver_Lineage",
            "Donor_Lineage",
        ],
        na_position="last",
    ).reset_index(drop=True)

    if output_file is not None:
        result_df.to_csv(folder / output_file, index=False)

    return result_df


A = None
B = None
C = None




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

for key, val, add in val_key:
    print(val)

    try:
        daft_folder_relations_to_newick(val)
    except Exception as e:
        print(f"{val} not found or failed")
        print(e)
        continue