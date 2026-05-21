# DAFT User Guide


DAFT, the **Discordant Attachment Frequency Tool**, is a Python package for detecting and characterizing introgression from a rooted species tree and a collection of gene trees.

This guide shows how to install DAFT, prepare input files, run the main analyses, and interpret the output.

---

## What DAFT Does

DAFT has two main components:

1. **DAFT Test**
   - Implemented in `DAFT_Test.py`
   - Counts discordant attachments between branches of a species tree across gene trees
   - Tests whether a lineage pair appears together more often than expected
   - Reports test counts, comparison counts, and Z-scores

2. **DAFT Direction**
   - Implemented in `DAFT_Direction.py`
   - Also referred to as `djiNNI`
   - Uses NNI-based reconciliation to infer the likely direction of introgression
   - Reports donor and receiver lineages
   - Produces an inferred introgression network

A typical DAFT workflow is:

```text
species tree + gene trees
        |
        v
DAFT_Test.py
        |
        v
significant lineage pairs
        |
        v
DAFT_Direction.py
        |
        v
donor/receiver inference + introgression network
```

---

## 1. Download DAFT

Clone the DAFT repository and move into the code directory:

```bash
git clone https://github.com/<your-username>/DAFT.git
cd DAFT
```

The main DAFT directory should contain:

```text
DAFT_Test.py
DAFT_Direction.py
DAFT_Transform.py
DAFT_produce_excel.py
DAFT_produce_excel_correction.py
clean_folder.py
excel_direction.py
DAFT_utils/
```

Run DAFT commands from this directory. The current version uses relative imports to access files in `DAFT_utils/`.

---

## 2. Install Dependencies

DAFT requires Python 3.

We recommend creating a clean environment:

```bash
conda create -n daft python=3.10
conda activate daft
```

Then install the required Python packages:

```bash
pip install numpy pandas ete3 openpyxl
```

The most important dependencies are:

```text
numpy
pandas
ete3
openpyxl
```

If you do not want to use conda, you can use a normal Python virtual environment:

```bash
python3 -m venv daft-env
source daft-env/bin/activate
pip install numpy pandas ete3 openpyxl
```

---

## 3. Prepare Input Files

DAFT needs two inputs:

1. A rooted species tree in Newick format
2. A CSV file containing gene trees in Newick format

---

## Species Tree

The species tree is passed directly with the `--sp` argument.

Example:

```text
(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));
```

In the command line, this becomes:

```bash
--sp "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));"
```

---

## Gene Tree CSV

The gene tree file is passed with the `--gt` argument.

Example:

```bash
--gt "test_data.csv"
```

The safest input format is a CSV file with one gene tree per row in the first column:

```csv
gt
"(A,((B,C),D));"
"((A,B),(C,D));"
```

The current code reads the CSV using `pandas.read_csv(...).to_numpy()` and uses the first column as the gene tree column, so gene trees should be placed in the first column.

---

## 4. Quick Start: Run the Full DAFT Workflow

The easiest way to run DAFT is to run `DAFT_Significance.py` and let it automatically call `DAFT_Direction.py`.

```bash
python3 DAFT_Significance.py \
  --sp "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));" \
  --gt "test_data.csv" \
  --output "output" \
  --sibling 1 \
  --excel 1 \
  --direction 1
```

This command:

- reads the species tree from `--sp`
- reads gene trees from `test_data.csv`
- uses `output` as the output label
- runs sibling comparisons with `--sibling 1`
- writes Excel output with `--excel 1`
- automatically runs direction inference with `--direction 1`

After the run finishes, DAFT organizes output files into:

```text
DAFT_results/
DAFT_extras/
```

---

## 5. Recommended First Run

For a first analysis, we recommend running the significance step first without automatic direction inference:

```bash
python3 DAFT_Significance.py \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --output "output" \
  --sibling 1 \
  --excel 1 \
  --direction 0
```

Then inspect:

```text
DAFT_results/DAFT_Significance_output.txt
DAFT_results/DAFT_Significance_output.xlsx
DAFT_results/branch_map.csv
```

After identifying significant lineage pairs, run direction inference manually:

```bash
python3 DAFT_Direction.py \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --lineages "[('B;', 'K;'), ('M;', 'E;')]" \
  --output "output"
```

This two-step workflow is the safest way to understand what DAFT is doing.

---

## 6. DAFT Significance Arguments

### `--sp`

Species tree in Newick format.

Example:

```bash
--sp "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));"
```

---

### `--gt`

CSV file containing gene trees.

Example:

```bash
--gt "test_data.csv"
```

---

### `--output`

Output label used in output file names.

Example:

```bash
--output "output"
```

This produces files such as:

```text
DAFT_Significance_output.txt
DAFT_Significance_output.xlsx
DAFT_Direction_output.txt
```

---

### `--sibling`

Controls whether DAFT includes sibling comparisons.

```text
1 = run sibling comparisons
0 = do not run sibling comparisons
```

Example:

```bash
--sibling 1
```

The sibling comparison gives an additional control comparison for the test lineage.

---

### `--excel`

Controls whether DAFT writes Excel output.

```text
1 = write Excel output
0 = do not write Excel output
```

Example:

```bash
--excel 1
```

When enabled, DAFT creates:

```text
DAFT_Significance_output.xlsx
```

---

### `--direction`

Controls whether DAFT automatically runs direction inference after the significance test.

```text
1 = automatically run DAFT_Direction.py
0 = only run DAFT_Significance.py
```

Example:

```bash
--direction 1
```

When this is enabled, DAFT selects significant lineage pairs and passes them to `DAFT_Direction.py`.

The automatic direction step currently uses:

```text
Z-score cutoff: -1.96
Minimum count cutoff: 6
```

This means DAFT Direction is only run on lineage pairs with both a strong enough Z-score and enough supporting gene-tree counts.

---

### `--correct`

Controls whether DAFT scales counts by the number of times each lineage or clade is observed in the gene tree set.

```text
1 = run corrected DAFT output
0 = run uncorrected DAFT output
```

Example:

```bash
--correct 1
```

Use `--correct 1` when some lineages or clades appear in fewer gene trees than others.

This can happen when:

- ILS is high
- an internal clade does not form in many gene trees
- some taxa are missing from some gene trees
- one comparison lineage is available many more times than another

Without correction, DAFT may compare raw attachment counts between lineages that were not equally available to be observed. This can make a test look significant simply because one comparison lineage was seen in many more gene trees than another.

With correction enabled, DAFT reports counts as:

```text
attachment_count / total_lineage_count
```

For example:

```text
5.00005/2000
```

means that the attachment was observed about 5 times, and the tested lineage was observed 2000 times across the input gene trees.

Another example:

```text
5e-05/779
```

means that the attachment count was 0, but DAFT added a very small pseudo-count so the corrected calculation can be performed safely.

When correction is enabled, Z-score columns may appear as:

```text
uncorrected_Z | corrected_Z
```

For example:

```text
-2.24|-0.73
```

means the raw count comparison looked significant, but after scaling by lineage availability, the corrected score was no longer significant.

---

### `--demography`

Optional demography file argument used in simulation or truth-checking workflows.

Example:

```bash
--demography "demography_M.txt"
```

For normal empirical analyses, this may not be needed.

---

## 7. Main Output Folders

DAFT organizes output into two folders.

---

## `DAFT_results/`

This folder contains the main output files.

Typical files include:

```text
DAFT_Significance_output.txt
DAFT_Significance_output.xlsx
DAFT_Direction_output.txt
branch_map.csv
important_output.csv
```

Some files only appear when the corresponding option is enabled.

For example:

- `DAFT_Significance_output.xlsx` appears when `--excel 1`
- `DAFT_Direction_output.txt` appears when `--direction 1` or when direction inference is run separately
- `branch_map.csv` maps DAFT branch labels back to the original species tree

---

## `DAFT_extras/`

This folder contains extra files created during the direction workflow.

These files are useful for inspecting gene-tree-level evidence and debugging the djiNNI step.

Examples may include files such as:

```text
djiNNI_out_B_K.csv
djiNNI_out_K_B.csv
```

The exact file names depend on the lineage pair and argument order.

---

## 8. Understanding `branch_map.csv`

DAFT labels branches internally. The file:

```text
DAFT_results/branch_map.csv
```

maps these internal branch labels back to the original species tree.

Example concept:

```text
Original branch: (B,U);
DAFT label:      U:39
```

Use this file whenever an output branch label is hard to interpret.

---

## 9. Understanding `DAFT_Significance_output.txt`

This is the main output file from `DAFT_Significance.py`.

It contains one section for each focal lineage.

A simplified example is:

```text
Species Tree = (((1,2),3),((4,5),6));
=====================================
Focal_lineage = ((4,5),6);

NNI_sp  Test_lineage  Test_count  comparison_uncle  uncle_count  Z-value-uncle  comparison_sibling  sibling_count  Z-value-sibling
====================================================================================================================================
2       1;            5           3;                0            -2.24          2;                  0              -2.24
2       (1,2);        1           -                 -            -              3;                  0              -1.00
------------------------------------------------------------------------------------------------------------------------------------
```

---

## Significance Output Columns

### `Focal_lineage`

The branch or clade being evaluated.

DAFT asks which other lineages attach discordantly to this focal lineage more often than expected.

---

### `Test_lineage`

The lineage being tested against the focal lineage.

---

### `NNI_sp`

The number of nearest-neighbor interchanges separating the focal lineage and test lineage in the species tree.

A larger value means the two lineages are farther apart in the species tree.

---

### `Test_count`

The number of gene trees in which the focal lineage and test lineage show the relevant discordant attachment.

In uncorrected output, this is a raw count.

In corrected output, this appears as:

```text
attachment_count / total_lineage_count
```

---

### `comparison_uncle`

The uncle, or avuncular lineage, used as a comparison for the test lineage.

The uncle comparison asks whether the test lineage attaches to the focal lineage more often than a closer comparison lineage does.

---

### `uncle_count`

The attachment count for the uncle lineage.

---

### `Z-value-uncle`

The Z-score comparing the test lineage count against the uncle lineage count.

More negative values indicate stronger evidence for introgression.

---

### `comparison_sibling`

The sibling lineage used as a second comparison.

This column is included when:

```bash
--sibling 1
```

---

### `sibling_count`

The attachment count for the sibling lineage.

---

### `Z-value-sibling`

The Z-score comparing the test lineage count against the sibling lineage count.

More negative values indicate stronger evidence for introgression.

---

## 10. Corrected Output Example

Correction is useful when one lineage is observed more often than another across the full gene tree set.

For example, suppose DAFT is comparing two lineages in a six-species tree:

```text
Species Tree = (((1,2),3),((4,5),6));
=====================================
Focal_lineage = ((4,5),6);

NNI_sp  Test_lineage  Test_count     comparison_uncle  uncle_count  Z-value-uncle  comparison_sibling  sibling_count  Z-value-sibling
=======================================================================================================================================
2       1;            5.00005/2000   3;                5e-05/215    -2.24|-0.73    2;                  5e-05/779      -2.24|-0.75
2       (1,2);        1.00005/803    -                 -            -|-            3;                  5e-05/779      -1.00|-0.98
---------------------------------------------------------------------------------------------------------------------------------------
```

The format:

```text
5.00005/2000
```

means the attachment was observed about 5 times, and that lineage was observed 2000 times in the full gene tree set.

The format:

```text
5e-05/779
```

means the attachment count was 0, but DAFT added a very small pseudo-count.

In the second row:

```text
Test_count:     1.00005/803
sibling_count:  5e-05/779
```

The test lineage was observed 803 times in the gene tree set, while its sibling comparison was observed 779 times.

Correction makes the comparison fair by scaling counts according to how often each lineage or clade was available to be observed.

The Z-score format:

```text
-2.24|-0.73
```

means:

```text
raw Z-score:       -2.24
corrected Z-score: -0.73
```

In this example, the raw test looks significant, but the corrected test does not. This means the apparent signal was partly caused by unequal lineage availability.

---

## 11. Interpreting Z-scores

DAFT uses negative Z-scores to identify unusually high discordant attachment frequencies.

A commonly used cutoff is:

```text
Z <= -1.96
```

More negative values indicate stronger evidence that the test lineage attaches to the focal lineage more often than expected.

In corrected output, focus on the corrected Z-score, which appears after the vertical bar:

```text
raw_Z | corrected_Z
```

Example:

```text
-2.24|-0.73
```

Here, the corrected value is:

```text
-0.73
```

---

## 12. Running DAFT Direction Separately

You can run `DAFT_Direction.py` directly when you already know which lineage pairs you want to test.

Example:

```bash
python3 DAFT_Direction.py \
  --sp "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));" \
  --gt "test_data.csv" \
  --lineages "[('B;', 'K;'), ('M;', 'E;')]" \
  --output "output"
```

This tests direction for:

```text
B; and K;
M; and E;
```

---

## 13. DAFT Direction Arguments

### `--sp`

Species tree in Newick format.

---

### `--gt`

CSV file containing gene trees.

---

### `--lineages`

A Python-style list of lineage pairs.

Example:

```bash
--lineages "[('B;', 'K;'), ('M;', 'E;')]"
```

Each tuple contains one pair of lineages to test.

---

### `--lineagesN`

Optional argument for bidirectional testing.

Example:

```bash
--lineagesN "[(['(1,2);', '3;'], '4;')]"
```

This provides extra sibling information used in the bidirectional test.

If this argument is not provided, the current code attempts to infer the required sibling information automatically.

---

### `--output`

Output label used in the output file name.

Example:

```bash
--output "output"
```

This produces:

```text
DAFT_Direction_output.txt
```

---

## 14. Understanding `DAFT_Direction_output.txt`

The direction output file usually contains these sections:

```text
SIGNIFICANT PAIRS
DATA TABLE1
DATA TABLE2(BIDIRECTIONAL)
INFERRED RELATIONS
NETWORK
```

---

## `SIGNIFICANT PAIRS`

This section lists the lineage pairs passed to DAFT Direction.

Example:

```text
BETWEEN 1; AND 4;           4 NNI APART
```

This means lineages `1;` and `4;` are 4 NNIs apart in the species tree.

---

## `DATA TABLE1`

This is the main direction inference table.

Important columns:

```text
Significant_Pairs
Total_gene_trees
Lineage1
Count1
Lineage2
Count2
Major_Moved
```

`Count1` and `Count2` summarize how often each lineage appears to move during the NNI transformations.

`Major_Moved` is the lineage with the higher movement count.

DAFT interprets `Major_Moved` as the likely receiver lineage.

---

## `DATA TABLE2(BIDIRECTIONAL)`

This section reports the bidirectional introgression test.

Important columns:

```text
Significant_Pairs
Minor_Moved
Minor_sibling
Total_gene_trees
Minor_sibling_count
Minor_moved_count
Z_score
```

This test checks whether there is support for introgression in the opposite direction.

A strongly negative Z-score supports a bidirectional signal.

---

## `INFERRED RELATIONS`

This section reports the inferred donor and receiver.

Example:

```text
BETWEEN 1; AND 4;        RECEIVER: 4; AND Donor: 1;         (BIDIRECTIONAL)
```

This means:

```text
Receiver = 4;
Donor    = 1;
```

If the line contains:

```text
(BIDIRECTIONAL)
```

then DAFT found evidence for introgression in both directions.

---

## `NETWORK`

The final section reports the inferred introgression network in Newick format.

Example:

```text
Network: (7,((6,(#H2,((#H5,5),((((4)#H1)#H3)#H6,#H4)))),(3,(#H1,((#H6,2),((((1)#H2)#H4)#H5,#H3))))));
```

This network can be visualized using a tree or network visualization tool that supports hybrid nodes or extended Newick-like notation.

---

## 15. Recommended Workflows

### Workflow A: significance only

Use this when you only want DAFT significance results:

```bash
python3 DAFT_Significance.py \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --output "output" \
  --sibling 1 \
  --excel 1 \
  --direction 0
```

Main outputs:

```text
DAFT_results/DAFT_Significance_output.txt
DAFT_results/DAFT_Significance_output.xlsx
DAFT_results/branch_map.csv
```

---

### Workflow B: corrected significance

Use this when lineages or clades are observed unequally across gene trees:

```bash
python3 DAFT_Significance.py \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --output "output" \
  --sibling 1 \
  --excel 1 \
  --correct 1 \
  --direction 0
```

Use the corrected Z-score when interpreting significance.

---

### Workflow C: manual direction inference

Use this after inspecting the significance output:

```bash
python3 DAFT_Direction.py \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --lineages "[('B;', 'K;'), ('M;', 'E;')]" \
  --output "output"
```

Main output:

```text
DAFT_Direction_output.txt
```

---

### Workflow D: full automatic workflow

Use this when you want DAFT to run significance and direction inference together:

```bash
python3 DAFT_Significance.py \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --output "output" \
  --sibling 1 \
  --excel 1 \
  --direction 1
```

Main outputs:

```text
DAFT_results/DAFT_Significance_output.txt
DAFT_results/DAFT_Significance_output.xlsx
DAFT_results/DAFT_Direction_output.txt
DAFT_results/branch_map.csv
DAFT_extras/
```

---

## 16. Troubleshooting

### `ModuleNotFoundError: No module named 'reconcILS'`

Run DAFT from the main DAFT directory:

```bash
cd DAFT
python3 DAFT_Significance.py ...
```

The current version uses relative paths to find `DAFT_utils/`.

---

### `ModuleNotFoundError: No module named 'ete3'`

Install `ete3`:

```bash
pip install ete3
```

---

### Excel output is missing

Make sure you used:

```bash
--excel 1
```

Also make sure `openpyxl` is installed:

```bash
pip install openpyxl
```

---

### Direction output is missing

Make sure you used:

```bash
--direction 1
```

If no pairs pass the significance and count filters, DAFT may not run direction inference automatically.

In that case, inspect the significance output and run `DAFT_Direction.py` manually with selected lineage pairs.

---

### Output files are not in the main folder

This is expected. DAFT moves outputs into:

```text
DAFT_results/
DAFT_extras/
```

Check those folders after the run.

---

### Branch names are hard to interpret

Use:

```text
DAFT_results/branch_map.csv
```

to map output labels back to the original species tree.

---

## 17. How We Recommend Reporting DAFT Results

When reporting DAFT results, include:

1. The species tree used.
2. The number of gene trees analyzed.
3. Whether `--correct 1` was used.
4. Whether sibling comparisons were used.
5. The Z-score threshold used for significance.
6. The significant lineage pairs.
7. The DAFT Direction donor/receiver inference, if direction inference was run.
8. Whether bidirectional introgression was inferred.
9. The final network, if reported.

Example language:

```text
We ran DAFT using a rooted species tree and 2,000 gene trees. 
Significance was assessed using avuncular and sibling comparisons. 
Lineage pairs with Z <= -1.96 and sufficient supporting counts were passed to DAFT Direction/djiNNI. 
Direction inference was performed using NNI reconciliation, and the lineage with the higher movement count was interpreted as the receiver.
```

---
