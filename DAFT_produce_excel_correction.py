import pandas as pd
import re
import argparse

def parse1():
    parser = argparse.ArgumentParser(description="DAFT_produce_excel_corrected")
    parser.add_argument('--output', type=str, help="Name of output file")
    args = parser.parse_args()
    return args

parser = parse1()
output = parser.output

input_file = "DAFT_Test_" + output + ".txt"
output_file = "DAFT_Test_" + output + ".xlsx"

columns = [
    "NNI_sp",
    "Test_lineage",
    "total_count",
    "comparison_uncle",
    "uncle_count",
    "Z_value_uncle",
    "comparison_sibling",
    "sibling_count",
    "Z_value_sibling",
]

with open(input_file) as f:
    lines = [l.rstrip("\n") for l in f]

blocks = []
current_focal = None
current_rows = []
col_starts = None 

for i, line in enumerate(lines):
    if line.startswith("Focal_lineage ="):
        if current_focal and current_rows:
            blocks.append((current_focal, current_rows))
        current_focal = line.replace("Focal_lineage =", "").strip()
        current_rows = []
        col_starts = None
        continue

    if line.startswith("=") or line.startswith("NNI_sp") or line.strip() == "":
        if line.startswith("NNI_sp"):
            col_starts = []
            last_pos = 0
            for col in columns:
                idx = line.find(col, last_pos)
                if idx == -1:
                    idx = last_pos
                col_starts.append(idx)
                last_pos = idx + len(col)
            col_starts.append(None) 
        continue

    if not re.match(r"^\d+", line):
        continue
    
    prt = []
    if col_starts:
        for j in range(len(col_starts) - 1):
            start = col_starts[j]
            end = col_starts[j+1]
            val = line[start:end].strip()
            prt.append(val)
        prt = prt[:9]
        while len(prt) < 9:
            prt.append("")
    else:

        prt = [p.strip() for p in re.split(r"\s{2,}", line)]
        while len(prt) < 9:
            prt.append("")
        prt = prt[:9]

    current_rows.append(prt)


if current_focal and current_rows:
    blocks.append((current_focal, current_rows))


with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    catalog = pd.DataFrame(
        [{"Focal_lineage": focal, "Comparison": f"{i}"} for i, (focal, row) in enumerate(blocks, 1)]
    )
    catalog.to_excel(writer, sheet_name="Catalog", index=False)

    for i, (focal, rows) in enumerate(blocks, 1):
        df = pd.DataFrame(rows, columns=columns)
        df.insert(0, "Focal_lineage", focal)
        df.to_excel(writer, sheet_name=f"Comparison_{i}", index=False)