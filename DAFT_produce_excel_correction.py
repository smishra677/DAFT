import pandas as pd
import re
import argparse


def parse1():
    parser = argparse.ArgumentParser(description="Convert DAFT txt output to Excel")
    parser.add_argument('--output', type=str, required=True, help="Name of output file")
    args = parser.parse_args()
    return args


def split_cols(line):
    """
    Split DAFT text columns.
    Handles both multiple spaces and tabs.
    """
    return [p.strip() for p in re.split(r"\t+|\s{2,}", line.strip()) if p.strip()]


args = parse1()
output = args.output

input_file = "DAFT_Test_" + output + ".txt"
output_file = "DAFT_Test_" + output + ".xlsx"

with open(input_file) as f:
    lines = [l.rstrip() for l in f]

blocks = []

current_focal = None
current_header = None
current_rows = []


def save_current_block():
    if current_focal and current_header and current_rows:
        blocks.append((current_focal, current_header, current_rows))


for line in lines:
    line = line.rstrip()

    if line.startswith("Focal_lineage ="):
        save_current_block()

        current_focal = line.replace("Focal_lineage =", "").strip()
        current_header = None
        current_rows = []
        continue

    # Detect selected columns from your DAFT text output
    if line.startswith("NNI_sp"):
        current_header = split_cols(line)
        continue

    # Skip separators and non-data lines
    if not current_header:
        continue

    if not line or line.startswith("=") or line.startswith("-"):
        continue

    # Data rows usually start with NNI_sp value
    if not re.match(r"^\d+", line.strip()):
        continue

    row = split_cols(line)

    # Pad missing columns
    while len(row) < len(current_header):
        row.append("")

    # Trim extra values if any
    row = row[:len(current_header)]

    current_rows.append(row)

save_current_block()

with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    catalog = pd.DataFrame(
        [
            {
                "Comparison": f"Comparison_{i}",
                "Focal_lineage": focal,
                "Rows": len(rows),
            }
            for i, (focal, header, rows) in enumerate(blocks, 1)
        ]
    )

    catalog.to_excel(writer, sheet_name="Catalog", index=False)

    for i, (focal, header, rows) in enumerate(blocks, 1):
        df = pd.DataFrame(rows, columns=header)
        df.insert(0, "Focal_lineage", focal)

        sheet_name = f"Comparison_{i}"
        df.to_excel(writer, sheet_name=sheet_name, index=False)

print(f"Excel file written to: {output_file}")