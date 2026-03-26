import pandas as pd
import re
import argparse


def parse1():
    parser = argparse.ArgumentParser(description="IQTree on Simphy and dupcoal")
   
    parser.add_argument('--output', type=str, help="Name of output file")
    
    args = parser.parse_args()
    return args

parser = parse1()
output =parser.output


input_file = "DAFT_Significance_"+output+".txt"
output_file = "DAFT_Significance_"+output+".xlsx"

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
    lines = [l.rstrip() for l in f]

blocks = []
current_focal = None
current_rows = []

for line in lines:
    if line.startswith("Focal_lineage ="):
        if current_focal and current_rows:
            blocks.append((current_focal, current_rows))
        current_focal = line.replace("Focal_lineage =", "").strip()
        current_rows = []
        continue

    if not re.match(r"^\d+", line):
        continue

    prt = [p.strip() for p in re.split(r"\s{2,}", line)]

    if len(prt) == 7:
        prt.insert(5, "")
        prt.append("")
    elif len(prt) == 8:
        if re.match(r"-?\d+(\.\d+)?$", prt[-1]):
            prt.insert(5, "")
        else:
            prt.append("")

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
