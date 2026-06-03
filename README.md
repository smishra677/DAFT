# DAFT  Version 1.0  User Guide


DAFT, the **Discordant Attachment Frequency Tool**, is a Python package for detecting and characterizing introgression from a rooted species tree and a collection of gene trees.

This guide shows how to install DAFT, prepare input files, run the main analyses, and interpret the output.

---

## What DAFT Does

DAFT has two main components:

1. **DAFT Test**
   - Implemented in `DAFT_Test.py`
   - Counts discordant attachments between branches of a species tree across gene trees
   - Tests whether a lineage pair appears together more often than expected
   - Reports test counts, comparison counts, and Z-scores for both the sister and avuncular tests

2. **DAFT Direction**
   - Implemented in `DAFT_Direction.py`
   - Also referred to as `djiNNI`
   - Uses NNI-based reconciliation to infer the likely direction of introgression
   - Reports donor and recipient lineages
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
donor/recipient inference + introgression network
```

---

## 1. Download DAFT

Clone the DAFT repository and move into the code directory:

```bash
git clone https://github.com/smishra677/DAFT.git
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
daft/
DAFT_utils/
```

Run DAFT commands from this directory. By default, DAFT looks for helper code in `./DAFT_utils`.

If you run DAFT from another location, pass the utility path explicitly using:

```bash
--path /path/to/DAFT_utils
```

---

## 2. Install DAFT

For full installation details, please read [`INSTALL.md`](INSTALL.md).

DAFT requires Python 3 and the following Python packages:

```text
pandas==2.2.2
numpy==1.26.4
matplotlib
igraph
ete3
openpyxl
rich
```

From the DAFT directory, install the dependencies:

```bash
python -m pip install -r requirements.txt
```

Then install DAFT:

```bash
python -m pip install .
```

After installation, the main DAFT commands are available from the command line:

```bash
daft-test --help
daft-direction --help
daft-transform --help
```

The installed commands automatically locate `DAFT_utils`.

DAFT can also be run directly from the original scripts:

```bash
python3 DAFT_Test.py
python3 DAFT_Direction.py
```

When running the original scripts from outside the DAFT directory, pass the utility path with:

```bash
--path /path/to/DAFT_utils
```

---

## 3. Prepare Input Files

DAFT needs two inputs:

1. A rooted species tree in Newick format
2. A gene tree file containing gene trees in Newick format

The species tree and gene trees should be rooted and bifurcating. Branch lengths are not required. Taxon names must be consistent between the species tree and the gene trees; gene trees may contain fewer taxa than the species tree, but they should not contain taxa absent from the species tree.

---

## Species Tree

The species tree can be passed directly with the `--sp` argument, or read from a file with the `--sp_file` argument.

Example:

```text
(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));
```

In the command line, this becomes:

```bash
--sp "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));"
```

The same species tree can also be saved in a plain text file and passed with:

```bash
--sp_file "species_tree.txt"
```

---

## Gene Tree File

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

DAFT also accepts a simple non-CSV tree-list file. In this format, each line contains one Newick gene tree:

```text
(A,((B,C),D));
((A,B),(C,D));
```

When a non-CSV file is passed with `--gt`, DAFT converts it internally to a CSV file with a single `gt` column.

---

## Input Validation and Runtime Checks

DAFT validates input files before running the main analysis. These checks are meant to make DAFT runs reproducible: a successful run should make clear which trees were analyzed, how they were rooted, which gene trees were excluded, and which output files should be inspected before interpreting results.

| Check | What DAFT expects | What to do if it fails |
|---|---|---|
| Species tree | One rooted, bifurcating Newick tree | Fix the species tree or provide it with `--sp_file` |
| Gene trees | One rooted, bifurcating Newick tree per row or line | Fix malformed Newick strings before running DAFT |
| Semicolons | Each complete tree ends with `;` | Add the final semicolon to incomplete tree strings |
| Taxon names | Gene-tree taxa are present in the species tree | Make taxon labels identical across files |
| Missing taxa | Allowed in gene trees and reported as warnings | Inspect the log and taxon report before interpretation |
| Duplicate taxa | Not allowed by default | Fix the tree, or intentionally skip duplicate-containing gene trees with `--ignore_duplication 1` |
| Rooting | Gene trees are rooted consistently with the species tree | Reroot with `--rooting OUTGROUP --forced 1`, or intentionally continue with `--allow_inconsistent_rooting 1` |
| Lineage-pair file | `DAFT_Direction.py` receives a valid Python-style list from `--lineages_file` | Fix brackets, commas, quotes, and tuple structure |

### Duplicate taxa

By default, DAFT stops if a tree contains duplicate taxa.

```text
(A,(B,B));
```

To skip gene trees with duplicate taxa instead of failing the full run, use:

```bash
--ignore_duplication 1
```

When duplicate-containing gene trees are skipped, DAFT writes a cleaned gene-tree file:

```text
<output>_duplicate_removed_gene_trees.csv
```

Use this option deliberately. Skipping trees changes the dataset used by the analysis and may change downstream significance and direction results.

### Rooting

DAFT expects rooted species trees and rooted gene trees. To report an outgroup without changing the trees, use:

```bash
--rooting A
```

This records the intended outgroup, but it does not reroot the input trees. To explicitly reroot the species tree and gene trees, use:

```bash
--rooting A --forced 1
```

Use this only if the listed taxa form the intended outgroup. If you intentionally want DAFT to continue without rerooting even when roots are inconsistent, use:

```bash
--allow_inconsistent_rooting 1
```

### Validation logs and audit files

DAFT writes validation and audit files using the output prefix supplied with `--output`.

| File | Produced by | Meaning |
|---|---|---|
| `<output>_DAFT_log.txt` | `DAFT_Test.py` | Run options, input paths, validation checks, warnings, and failures |
| `<output>_djiNNI_log.txt` | `DAFT_Direction.py` | Direction-stage validation and djiNNI run log |
| `<output>_taxon_report.csv` | validation stage | Taxon presence across species tree and gene trees |
| `<output>_rooted_gene_trees.csv` | when `--forced 1` is used | Gene trees after forced rerooting |
| `<output>_duplicate_removed_gene_trees.csv` | when `--ignore_duplication 1` removes trees | Gene trees retained after duplicate-containing trees are skipped |

After cleanup, these validation files are placed in:

```text
DAFT_log_<output>/
```

The taxon report contains:

```text
original_taxon
in_species_tree
gene_tree_count
daft_internal_name
```

Use the taxon report to confirm that all taxon names were interpreted as intended.

---

## 4. Quick Start: Run the Full DAFT Workflow

The main DAFT workflow starts with `daft-test`. When `--direction 1` is used, DAFT also runs `daft-direction` on the significant lineage pairs.

```bash
daft-test \
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

The same workflow can also be run with the original script:

```bash
python3 DAFT_Test.py \
  --sp "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));" \
  --gt "test_data.csv" \
  --path "./DAFT_utils" \
  --output "output" \
  --sibling 1 \
  --excel 1 \
  --direction 1
```

After the run finishes, DAFT organizes output files into folders named with the output label:

```text
DAFT_results_<output>/
DAFT_extras_<output>/
DAFT_log_<output>/
```

---

## 5. Recommended First Run

For a first analysis, run the significance step first without automatic direction inference:

```bash
daft-test \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --output "output" \
  --sibling 1 \
  --excel 1 \
  --direction 0
```

Then inspect:

```text
DAFT_results_<output>/DAFT_Test_<output>.txt
DAFT_results_<output>/DAFT_Test_<output>.xlsx
DAFT_results_<output>/branch_map.csv
```

After identifying significant lineage pairs, run direction inference manually:

```bash
daft-direction \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --lineages_file "lineages.txt" \
  --verbose 1 \
  --output "output"
```

This two-step workflow is the safest way to understand what DAFT is doing.

The original script commands are also supported:

```bash
python3 DAFT_Test.py \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --path "./DAFT_utils" \
  --output "output" \
  --sibling 1 \
  --excel 1 \
  --direction 0

python3 DAFT_Direction.py \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --path "./DAFT_utils" \
  --lineages_file "lineages.txt" \
  --verbose 1 \
  --output "output"
```

---

## 6. DAFT Test Arguments

Use this table as the argument reference for `daft-test` or `python3 DAFT_Test.py`.

| Argument | Required? | Default | Input format / allowed values | Example | Meaning / output effect |
|---|---:|---:|---|---|---|
| `--sp` | Required if `--sp_file` is not used | `None` | Species tree as a Newick string | `--sp "(A,(B,C));"` | Supplies the species tree directly on the command line. DAFT adds a final semicolon if needed. |
| `--sp_file` | Required if `--sp` is not used | `None` | Path to a text file containing one species-tree Newick string | `--sp_file "species_tree.txt"` | Reads the species tree from a file. Use either `--sp` or `--sp_file`. |
| `--gt` | Yes | none | Path to a CSV file or a simple tree-list file | `--gt "test_data.csv"` | Supplies the gene tree set. For CSV input, DAFT expects a gene-tree column; for a simple tree-list file, DAFT converts it to CSV internally. |
| `--path` | No | `./DAFT_utils` | Path to the DAFT utility folder | `--path "./DAFT_utils"` | Lets DAFT find helper modules. Change this when running from outside the repository root or when `DAFT_utils/` is stored elsewhere. |
| `--output` | No, but strongly recommended | `None` | Output label string | `--output "output"` | Used in output names such as `DAFT_Test_<output>.txt`, `DAFT_Test_<output>.xlsx`, and `DAFT_Direction_<output>.txt`. |
| `--sibling` | No | `0` | `1` or `0` | `--sibling 1` | Controls whether DAFT includes sibling comparisons. `1` runs sibling comparisons; `0` skips them. |
| `--excel` | No | `0` | `1` or `0` | `--excel 1` | Controls whether DAFT writes `DAFT_Test_<output>.xlsx`. The Excel output is useful when the text output is too large to inspect comfortably. |
| `--direction` | No | `0` | `1` or `0` | `--direction 1` | Controls whether DAFT automatically runs `DAFT_Direction.py` after DAFT Test. `1` runs direction inference on significant pairs; `0` only runs DAFT Test. |
| `--correct` | No | `0` | `1` or `0` | `--correct 1` | Controls whether DAFT reports corrected output. `1` adds appearance columns and corrected Z-score columns; `0` reports uncorrected output only. |
| `--demography` | No | `None` | Path to a demography file | `--demography "demography.txt"` | Optional parser argument for demography-related workflows. Standard DAFT Test runs do not require it. |
| `--rooting` | No | `None` | Outgroup taxon, or comma-separated outgroup taxa | `--rooting A` | Records the intended outgroup for rooting checks. When used with `--forced 1`, DAFT reroots the species tree and gene trees using this outgroup. |
| `--forced` | No | `0` | `1` or `0` | `--rooting A --forced 1` | Controls explicit rerooting. `1` reroots using `--rooting`; `0` does not reroot. |
| `--allow_inconsistent_rooting` | No | `0` | `1` or `0` | `--allow_inconsistent_rooting 1` | Controls whether DAFT continues when gene-tree roots are inconsistent with the species-tree root. Use `1` only when the inconsistency is intentional. |
| `--ignore_duplication` | No | `0` | `1` or `0` | `--ignore_duplication 1` | Controls duplicate-taxon handling in gene trees. `1` skips duplicate-containing gene trees and writes `<output>_duplicate_removed_gene_trees.csv`; `0` stops when duplicates are detected. |
| `--random_seed` | No | `42` | Integer | `--random_seed 42` | Sets the random seed used in the djiNNI direction step when tied movement counts must be resolved. |
| `--verbose` | No | `0` | `1` or `0` | `--verbose 1` | Controls validation and progress messages. `1` prints more information; `0` is quieter. |

The automatic direction step uses the significance filters currently implemented in DAFT: a Z-score cutoff of `-1.96` and a minimum count cutoff of `6`. This means DAFT Direction is only run automatically on lineage pairs with a strong enough Z-score and enough supporting gene-tree counts.

---

## 7. Main Output Folders

DAFT organizes output into folders named with the output label. The exact files present depend on the arguments used.

| Folder | File or file pattern | Appears when | Purpose |
|---|---|---|---|
| `DAFT_results_<output>/` | `DAFT_Test_<output>.txt` | Standard DAFT Test run | Main fixed-width text output from `DAFT_Test.py`. |
| `DAFT_results_<output>/` | `DAFT_Test_<output>.xlsx` | `--excel 1` | Excel version of the DAFT Test output, useful for large species trees. |
| `DAFT_results_<output>/` | `DAFT_Direction_<output>.txt` | `--direction 1`, or a separate DAFT Direction run | Main fixed-width text output from `DAFT_Direction.py`. |
| `DAFT_results_<output>/` | `branch_map.csv` | Standard DAFT run | Maps DAFT numeric branch IDs back to species-tree branches. Use this file to interpret numeric lineage IDs in the output. |
| `DAFT_extras_<output>/` | `djiNNI_<output>_<lineage1_id>_<lineage2_id>.csv` | Direction inference | Intermediate djiNNI comparison file. Numeric IDs are used to avoid very long file names when lineage descriptions contain full species names. |
| `DAFT_extras_<output>/` | Other `introgression*`, `rev*`, `result*`, `list_*`, or `Summary*` CSV files | Direction inference / djiNNI internals | Extra files useful for inspecting gene-tree-level evidence and debugging the djiNNI step. |
| `DAFT_log_<output>/` | `<output>_DAFT_log.txt` | Validation / DAFT Test | Human-readable validation and run log. |
| `DAFT_log_<output>/` | `<output>_djiNNI_log.txt` | Direction inference | Human-readable djiNNI progress log. |
| `DAFT_log_<output>/` | `<output>_taxon_report.csv` | Validation | CSV report showing taxon presence across the species tree and gene trees. |
| `DAFT_log_<output>/` | `<output>_rooted_gene_trees.csv` | `--rooting OUTGROUP --forced 1` | Gene trees after explicit rerooting. |
| `DAFT_log_<output>/` | `<output>_duplicate_removed_gene_trees.csv` | `--ignore_duplication 1` removes one or more trees | Cleaned gene-tree file after duplicate-containing gene trees are skipped. |

Because DAFT now reports full species names in lineage descriptions, using full lineage strings directly in intermediate djiNNI file names can make file names too long. To avoid this, DAFT uses the same numeric branch IDs used in the main DAFT output and in `branch_map.csv` when naming extra djiNNI files. For example, if the output prefix is `out`, the intermediate comparison between lineages `2` and `3` is written as `djiNNI_out_2_3.csv` rather than a file name containing the full Newick or full species-name description for both lineages.

---

## 8. Understanding `branch_map.csv`

DAFT labels branches internally. The file:

```text
DAFT_results_<output>/branch_map.csv
```

maps these internal branch labels back to the original species tree. Use this file whenever an output branch label is hard to interpret.

### `branch_map.csv` columns

| Column | Meaning |
|---|---|
| `From` | Parent branch or clade in the labelled species tree. For the root, this can be the same as `To`. |
| `To` | Child branch, terminal taxon, or clade assigned a DAFT numeric ID. |
| `id` | Numeric DAFT branch ID used in the main output tables, DAFT Direction tables, network labels, and extra djiNNI file names. |

Example row interpretation:

| `From` | `To` | `id` | Meaning |
|---|---|---:|---|
| `(B,U);` | `U;` | `39` | Branch ID `39` refers to the branch leading to `U` from the `(B,U)` clade. |

---

### `<output>_taxon_report.csv` columns

This validation file is written to `DAFT_log_<output>/` and reports how taxa were seen by DAFT.

| Column | Meaning |
|---|---|
| `original_taxon` | Taxon name as read from the species tree and gene trees. |
| `in_species_tree` | `True` if the taxon appears in the species tree; `False` otherwise. |
| `gene_tree_count` | Number of gene trees containing the taxon. |
| `daft_internal_name` | Internal DAFT name assigned to the taxon. This is usually the same as the original taxon name. |

---

### Rooted and duplicate-filtered gene-tree CSV columns

Files such as `<output>_rooted_gene_trees.csv` and `<output>_duplicate_removed_gene_trees.csv` use this simple format:

| Column | Meaning |
|---|---|
| `gt` | Gene tree Newick string after the relevant validation step. |

---

## 9. Understanding `DAFT_Test_output.txt`

This is the main output file from `DAFT_Test.py`.

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

## Test Output Format

`DAFT_Test_<output>.txt` is organized as one block per focal lineage. Each block starts with a `Focal_lineage = ...` line, followed by a fixed-width table.

### DAFT Test block structure

| Output item | Appears when | Meaning |
|---|---|---|
| `Species Tree = ...` | Once at the top of the file | The labelled species tree used by DAFT. Branch labels in this tree correspond to the numeric IDs used in the output tables. |
| `Focal_lineage = ...` | At the start of each focal-lineage block | The branch or clade currently being evaluated. DAFT asks which other lineages attach discordantly to this focal lineage more often than expected. |
| Header row beginning with `NNI_sp` | In each focal-lineage block | Column names for the fixed-width output table. The exact columns depend on whether correction and sibling comparisons are enabled. |
| Data rows | In each focal-lineage block | One row per test lineage compared against the current focal lineage. |
| Separator row of `---` characters | After each focal-lineage block | Marks the end of that focal-lineage block. |

### DAFT Test table columns

| Column | Appears when | Meaning |
|---|---|---|
| `NNI_sp` | Always | Number of nearest-neighbor interchanges separating the focal lineage and test lineage in the species tree. Larger values mean the lineages are farther apart in the species tree. |
| `Test_lineage` | Always | The lineage being tested against the current focal lineage. |
| `Test_appearance` | `--correct 1` | Number of gene trees in which the test lineage or clade was available to be observed. |
| `Test_count` | Always | Observed number of gene trees in which the focal lineage and test lineage show the relevant discordant attachment. This remains an observed count, not a corrected count. |
| `comparison_uncle` | Always | The uncle, or avuncular, lineage used as a comparison for the test lineage. A blank or `-` value means no valid uncle comparison was available for that row. |
| `uncle_appearance` | `--correct 1` | Number of gene trees in which the uncle lineage or clade was available to be observed. |
| `uncle_count` | Always | Observed attachment count for the uncle comparison lineage. |
| `Z-value-uncle` | Always | Raw Z-score comparing the test lineage count against the uncle lineage count. More negative values indicate stronger evidence for introgression. |
| `Z-value-uncle_corrected` | `--correct 1` | Corrected uncle-comparison Z-score after scaling by lineage availability. |
| `comparison_sibling` | `--sibling 1` | Sibling lineage used as a second comparison. A blank or `-` value means no valid sibling comparison was available for that row. |
| `sibling_count` | `--sibling 1` | Observed attachment count for the sibling comparison lineage. |
| `sibling_appearance` | `--correct 1` and `--sibling 1` | Number of gene trees in which the sibling lineage or clade was available to be observed. |
| `Z-value-sibling` | `--sibling 1` | Raw Z-score comparing the test lineage count against the sibling lineage count. More negative values indicate stronger evidence for introgression. |
| `Z-value-sibling_corrected` | `--correct 1` and `--sibling 1` | Corrected sibling-comparison Z-score after scaling by lineage availability. |

In corrected output, DAFT keeps observed counts and availability counts in separate columns. For example, `Test_count` is the observed attachment count, while `Test_appearance` is the number of gene trees in which that test lineage or clade was available to be observed.

---

## 10. Corrected Output Example

Correction is useful when one lineage is observed more often than another across the full gene tree set. Corrected output has the same focal-lineage block structure as uncorrected output, but it adds appearance columns and corrected Z-score columns.

Example corrected DAFT Test block:

```text
Focal_lineage = 11
NNI_sp  Test_lineage  Test_appearance  Test_count       comparison_uncle  uncle_appearance  uncle_count  Z-value-uncle       Z-value-uncle_corrected  comparison_sibling  sibling_count  sibling_appearance  Z-value-sibling  Z-value-sibling_corrected
==================================================================================================================================================================================================================================================
0       6             1584             1264                              0                 0            -                   -                                            0              0                   -                -
1       8             1613             156                               0                 0            -                   -                        7                   138            2000                -1.05            -2.73
1       7             2000             138                               0                 0            -                   -                        8                   156            1613                1.05             2.73
1       12            2000             163                               0                 0            -                   -                                            0              0                   -                -
2       13            2000             19                12                2000              163          10.67               10.67                                        0              0                   -                -
2       10            2000             33                7                 2000              138          8.03                8.03                     9                   31             2000                -0.25            -0.25
2       9             2000             31                7                 2000              138          8.23                8.23                     10                  33             2000                0.25             0.25
3       14            1524             2                 13                2000              19           3.71                3.07                                         0              0                   -                -
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```

In this corrected output, `Test_count`, `uncle_count`, and `sibling_count` are still the observed attachment counts.

The appearance columns report how often each lineage or clade was available to be observed:

| Appearance column | Count column it helps interpret | Meaning |
|---|---|---|
| `Test_appearance` | `Test_count` | Number of gene trees in which the test lineage or clade was available to be observed. |
| `uncle_appearance` | `uncle_count` | Number of gene trees in which the uncle comparison lineage or clade was available to be observed. |
| `sibling_appearance` | `sibling_count` | Number of gene trees in which the sibling comparison lineage or clade was available to be observed. |

For example, in the row with `Test_lineage = 8`, the test lineage was available in `1613` gene trees and attached to the focal lineage in `156` gene trees. The sibling lineage `7` was available in `2000` gene trees and had `138` observed attachments. The raw sibling Z-score is `-1.05`, while the corrected sibling Z-score is `-2.73`.

The raw and corrected Z-scores are printed in separate columns:

| Raw Z-score column | Corrected Z-score column |
|---|---|
| `Z-value-uncle` | `Z-value-uncle_corrected` |
| `Z-value-sibling` | `Z-value-sibling_corrected` |

Correction makes the comparison fairer by scaling counts according to how often each lineage or clade was available to be observed.

***IMPORTANT NOTE: biologically significant results should have significant Z-scores BOTH before and after correction. We do not recommend using results that are only significant after correction.***

---

## 11. Interpreting Z-scores

DAFT uses negative Z-scores to identify unusually high discordant attachment frequencies.

A commonly used cutoff (for a nominal p-value < 0.05) is:

```text
Z <= -1.96
```

More negative values indicate stronger evidence that the test lineage attaches to the focal lineage more often than expected.

In corrected output, DAFT also reports corrected Z-scores in separate columns rather than combining raw and corrected values into one field.

For example:

```text
Z-value-uncle            -2.24
Z-value-uncle_corrected  -0.73
```

Here, the corrected value is `-0.73`.

---

## 12. Running DAFT Direction Separately

Run DAFT Direction when you already know which lineage pairs you want to test.

Create a plain text file containing the lineage pairs:

```text
[('B;', 'K;'), ('M;', 'E;')]
```

Then pass this file with `--lineages_file`:

```bash
daft-direction \
  --sp "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));" \
  --gt "test_data.csv" \
  --lineages_file "lineages.txt" \
  --verbose 1 \
  --output "output"
```

The same command can also be run with the original script:

```bash
python3 DAFT_Direction.py \
  --sp "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));" \
  --gt "test_data.csv" \
  --path "./DAFT_utils" \
  --lineages_file "lineages.txt" \
  --verbose 1 \
  --output "output"
```

This tests direction for:

```text
B; and K;
M; and E;
```

---

## 13. DAFT Direction Arguments

Use this table as the argument reference for `daft-direction` or `python3 DAFT_Direction.py`.

| Argument | Required? | Default | Input format / allowed values | Example | Meaning / output effect |
|---|---:|---:|---|---|---|
| `--sp` | Required if `--sp_file` is not used | `None` | Species tree as a Newick string | `--sp "(A,(B,C));"` | Supplies the species tree directly on the command line. |
| `--sp_file` | Required if `--sp` is not used | `None` | Path to a text file containing one species-tree Newick string | `--sp_file "species_tree.txt"` | Reads the species tree from a file. Use either `--sp` or `--sp_file`. |
| `--gt` | Yes for normal runs | none | Path to a CSV file or simple tree-list file | `--gt "test_data.csv"` | Supplies the gene tree set used for djiNNI direction inference. |
| `--path` | No | `./DAFT_utils` | Path to the DAFT utility folder | `--path "./DAFT_utils"` | Lets DAFT Direction find helper modules. Change this when running from outside the repository root or when `DAFT_utils/` is stored elsewhere. |
| `--verbose` | No | `1` | `1` or `0` | `--verbose 1` | Controls progress messages in DAFT Direction. `1` prints more information; `0` is quieter. |
| `--lineages_file` | Yes | none | Path to a text file containing a Python-style list of lineage tuples | `--lineages_file "lineages.txt"` | Supplies the significant lineage pairs to test. Example file contents: `[('B;', 'K;'), ('M;', 'E;')]`. |
| `--lineagesN_file` | No | `None` | Path to a text file containing non-unique or bidirectional lineage tuples | `--lineagesN_file "lineagesN.txt"` | Optional bidirectional-lineage input. Manual DAFT Direction runs usually do not need this file; automatic direction runs from DAFT Test create and pass this when needed. |
| `--output` | No, but strongly recommended | `None` | Output label string | `--output "output"` | Used in output names such as `DAFT_Direction_<output>.txt` and intermediate djiNNI files. |
| `--rooting` | No | `None` | Outgroup taxon, or comma-separated outgroup taxa | `--rooting A` | Records the intended outgroup for rooting checks. When used with `--forced 1`, DAFT reroots the species tree and gene trees using this outgroup. |
| `--forced` | No | `0` | `1` or `0` | `--rooting A --forced 1` | Controls explicit rerooting. `1` reroots using `--rooting`; `0` does not reroot. |
| `--allow_inconsistent_rooting` | No | `0` | `1` or `0` | `--allow_inconsistent_rooting 1` | Controls whether DAFT continues when gene-tree roots are inconsistent with the species-tree root. Use `1` only when the inconsistency is intentional. |
| `--cache_hash` | No | `None` | Hash string | `--cache_hash "<hash>"` | Input hash used to separate djiNNI cache results by dataset. DAFT Test passes this automatically during automatic direction inference; manual users usually do not need to set it. |
| `--ignore_duplication` | No | `0` | `1` or `0` | `--ignore_duplication 1` | Controls duplicate-taxon handling in gene trees. `1` skips duplicate-containing gene trees; `0` stops when duplicates are detected. |
| `--random_seed` | No | `42` | Integer | `--random_seed 42` | Sets the random seed used by djiNNI when tied movement counts must be resolved. |

---

## djiNNI Cache

`DAFT_Direction.py` / `djiNNI` may create a cache folder in the main DAFT working directory:

`djiNNI_cache/`

This folder stores previously computed reconciliation results for gene trees. The purpose of the cache is to make repeated `djiNNI` runs faster by reusing results that were already calculated.

The cache is tied to the gene tree order used during the run. This means the cache is only safe to reuse when you are running DAFT on the same gene tree dataset and the gene trees are in the same order as the previous run.

### When to keep `djiNNI_cache/`

You can keep the cache folder if you are rerunning DAFT on the same gene tree file and the gene tree order has not changed.

For example, it is safe to keep the cache if:

- `test_data.csv` is the same file
- gene trees are in the same order
- only the lineage pair or output name changed

In this case, keeping `djiNNI_cache/` can save time because DAFT can reuse previously computed reconciliation results.

### When to delete or move `djiNNI_cache/`

Delete or move the cache folder before running DAFT if you are using a different dataset.

You should also delete or move the cache if:

- you changed the gene tree CSV file
- you added gene trees
- you removed gene trees
- you reordered gene trees
- you are running a completely different analysis
- you are not sure whether the current cache matches the current gene tree file

If the cache folder was created from a different gene tree file or a different gene tree order, DAFT may reuse reconciliation results that do not belong to the current dataset.

---

## 14. Understanding `DAFT_Direction_output.txt`

## DAFT Direction
```text
SIGNIFICANT PAIRS
================================================================================
BETWEEN 26 AND 23                           4 NNI APART 
BETWEEN 26 AND 22                           4 NNI APART 
BETWEEN 38 AND 42                           1 NNI APART 
BETWEEN 17 AND 25                           1 NNI APART 
********************************************************************************



DATA TABLE1 (ONLY NNI > 1) 
====================================================================================
Significant_Pairs  Total_gene_trees  Lineage1  Count1  Lineage2  Count2  Major_Moved
====================================================================================
26 AND 23          35                26        35      23        0       26         
26 AND 22          11                26        11      22        0       26         
************************************************************************************



DATA TABLE2 (BIDIRECTIONAL) (ONLY NNI > 1) 
================================================================================================================
Significant_Pairs  Minor_Moved  Minor_sibling  Total_gene_trees  Minor_sibling_count  Minor_moved_count  Z_score           
================================================================================================================
26 AND 23          23           22             11                0                    0                  -                 
26 AND 22          22           23             35                0                    0                  -                 
****************************************************************************************************************



INFERRED RELATIONS (ONLY NNI > 1) :
================================================================================
26 AND 23                 RECEIVER:26 AND Donor:23
26 AND 22                 RECEIVER:26 AND Donor:22
********************************************************************************


================================================================================
Network (ONLY NNI > 1) :  (((((((P[&label=7],(R[&label=9],Q[&label=10])),O[&label=11]),S[&label=12]),T[&label=13]),(((((F[&label=19],G[&label=20]),((#H2,I[&label=22]),(#H1,C[&label=23]))),((((H[&label=26])#H1)#H2,D[&label=27]),E[&label=28])),J[&label=29]),((N[&label=32],M[&label=33]),(L[&label=35],K[&label=36])))),((U[&label=39],B[&label=40]),(V[&label=42],W[&label=43]))),A[&label=44]);
********************************************************************************
```



The direction output file usually contains these sections:

| Output section | Appears when | Purpose |
|---|---|---|
| `SIGNIFICANT PAIRS` | Always | Lists the significant lineage pairs passed to DAFT Direction and the NNI distance between each pair. |
| `DATA TABLE1 (ONLY NNI > 1)` | When at least one significant pair is more than 1 NNI apart | Summarizes which lineage moved more often during NNI transformations. This table is used to infer the likely receiver. |
| `DATA TABLE2 (BIDIRECTIONAL) (ONLY NNI > 1)` | When at least one significant pair is more than 1 NNI apart | Tests whether there is evidence for introgression in the opposite direction. |
| `INFERRED RELATIONS (ONLY NNI > 1)` | When at least one `NNI > 1` relation is inferred | Reports receiver and donor lineages inferred from the direction tables. |
| `Network (ONLY NNI > 1)` | When a network can be written from inferred `NNI > 1` relations | Reports the inferred introgression network in Newick-like format. |

---

## `SIGNIFICANT PAIRS`

This section lists the lineage pairs passed to DAFT Direction and the NNI distance between the two lineages in the species tree.

Although the output is printed as sentences, each line can be read as a two-column table:

| Significant pair | NNI distance |
|---|---:|
| `26 AND 23` | `4 NNI APART` |
| `26 AND 22` | `4 NNI APART` |
| `38 AND 42` | `1 NNI APART` |
| `17 AND 25` | `1 NNI APART` |

Using the example above, this line:

```text
BETWEEN 26 AND 23                           4 NNI APART
```

means that lineages `26` and `23` passed the DAFT significance filters and are 4 NNI moves apart in the species tree.

The `SIGNIFICANT PAIRS` section can include both `NNI > 1` pairs and `1 NNI APART` pairs. In the example above, `38 AND 42` and `17 AND 25` are listed because they passed the significance filters, but they are exactly 1 NNI apart.

---

## `DATA TABLE1 (ONLY NNI > 1)`

This is the main direction inference table for significant pairs that are more than 1 NNI apart.

| Column | Meaning |
|---|---|
| `Significant_Pairs` | Pair of lineage IDs being tested for direction. |
| `Total_gene_trees` | Number of gene trees contributing to this direction comparison. |
| `Lineage1` | First lineage in the significant pair. |
| `Count1` | Number of times `Lineage1` was counted as moving during the NNI transformations. |
| `Lineage2` | Second lineage in the significant pair. |
| `Count2` | Number of times `Lineage2` was counted as moving during the NNI transformations. |
| `Major_Moved` | Lineage with the higher movement count. DAFT interprets this as the likely receiver lineage. |

Using the example above:

```text
26 AND 23          35                26        35      23        0       26
```

means that 35 gene trees supported the comparison between lineages `26` and `23`. Lineage `26` moved 35 times and lineage `23` moved 0 times, so DAFT identifies `26` as the major moved lineage and therefore the likely receiver.

Pairs that are exactly 1 NNI apart are not shown in this table because this table is labelled `ONLY NNI > 1`.

---

## `DATA TABLE2 (BIDIRECTIONAL) (ONLY NNI > 1)`

This section reports the bidirectional introgression test for the `NNI > 1` pairs.

| Column | Meaning |
|---|---|
| `Significant_Pairs` | Pair of lineage IDs being tested. |
| `Minor_Moved` | The lineage with the lower movement count in `DATA TABLE1`. This is the inferred donor in the primary direction test. |
| `Minor_sibling` | Sibling/comparison lineage used to test whether the opposite direction also has support. |
| `Total_gene_trees` | Number of gene trees contributing to the bidirectional comparison. |
| `Minor_sibling_count` | Attachment or movement count for the minor sibling comparison lineage. |
| `Minor_moved_count` | Attachment or movement count for the minor moved lineage. |
| `Z_score` | Z-score for the bidirectional test. A strongly negative value supports a bidirectional signal; `-` means the value was not available or not applicable. |

Using the example above:

```text
26 AND 23          23           22             11                0                    0                  -
```

means that DAFT checked whether there was evidence for the opposite direction involving minor moved lineage `23` and its comparison sibling `22`. The `Z_score` is `-`, so this row does not provide evidence for bidirectional introgression.

---

## `INFERRED RELATIONS (ONLY NNI > 1)`

This section reports the inferred donor and receiver for each `NNI > 1` pair.

Each line can be read as:

| Output field | Meaning |
|---|---|
| Pair label, such as `26 AND 23` | The significant pair being interpreted. |
| `RECEIVER:<id>` | The inferred receiver lineage. This usually matches `Major_Moved` from `DATA TABLE1`. |
| `Donor:<id>` | The inferred donor lineage. This is the other lineage in the pair. |
| `(BIDIRECTIONAL)` | Optional label printed when the bidirectional test also supports the opposite direction. |

Using the example above:

```text
26 AND 23                 RECEIVER:26 AND Donor:23
```

means:

```text
Receiver = 26
Donor    = 23
```

This follows from `DATA TABLE1`, where lineage `26` had the higher movement count and was identified as `Major_Moved`.

---

## `Network (ONLY NNI > 1)`

The final section reports the inferred introgression network in Newick-like format.

| Output field | Meaning |
|---|---|
| `Network (ONLY NNI > 1)` | Indicates that the network is built only from significant pairs more than 1 NNI apart. |
| Newick-like network string | The inferred introgression network, including hybrid-node labels when direction relations are inserted. |
| `[&label=<id>]` annotations | Branch or taxon labels corresponding to DAFT's numeric lineage IDs. Use `branch_map.csv` to map these labels back to the species-tree branches. |

In the example above, the network is built from the inferred `NNI > 1` relations:

```text
26 AND 23                 RECEIVER:26 AND Donor:23
26 AND 22                 RECEIVER:26 AND Donor:22
```

This is why the network section is labelled `ONLY NNI > 1`. The `1 NNI APART` pairs from `SIGNIFICANT PAIRS`, such as `38 AND 42` and `17 AND 25` in the example above, are reported as significant pairs but are not included in `DATA TABLE1`, `DATA TABLE2`, `INFERRED RELATIONS`, or the final network.

The network can be visualized using a tree or network visualization tool that supports hybrid nodes or extended Newick-like notation.

---

## 15. Recommended Workflows

### Workflow A: significance only

Use this when you only want DAFT significance results:

```bash
daft-test \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --output "output" \
  --sibling 1 \
  --excel 1 \
  --direction 0
```

Main outputs:

```text
DAFT_results_<output>/DAFT_Test_<output>.txt
DAFT_results_<output>/DAFT_Test_<output>.xlsx
DAFT_results_<output>/branch_map.csv
```

---

### Workflow B: corrected significance

After running DAFT in normal mode, users can check to see if results are still significant after correction. This mode corrects for cases in which branches or clades are observed are missing:

```bash
daft-test \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --output "output" \
  --sibling 1 \
  --excel 1 \
  --correct 1 \
  --direction 0
```

Use the uncorrected and corrected Z-scores together when interpreting significance.

---

### Workflow C: manual direction inference

Use this after inspecting the significance output:

```bash
daft-direction \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --lineages_file "lineages.txt" \
  --verbose 1 \
  --output "output"
```

Main output:

```text
DAFT_Direction_<output>.txt
```

---

### Workflow D: full automatic workflow

Use this when you want DAFT to run significance and direction inference together:

```bash
daft-test \
  --sp "<species_tree_newick>" \
  --gt "test_data.csv" \
  --output "output" \
  --sibling 1 \
  --excel 1 \
  --direction 1
```

The original script commands are still supported. Add `--path "./DAFT_utils"` when using `python3 DAFT_Test.py` or `python3 DAFT_Direction.py` directly.

---

## 16. Rich Progress Dashboard and Verbose Output

`DAFT_Direction.py` and `DAFT_Transform.py` can print a rich terminal progress display during DAFT runs when `--verbose 1` is used. This display is meant to help users monitor long analyses while gene trees are being reconciled, transformed, or loaded from cache.

`DAFT_Transform.py` uses `rich` to show a terminal progress dashboard during the djiNNI transformation step.

<img width="976" height="527" alt="image" src="https://github.com/user-attachments/assets/8f8244ed-3754-4e08-89b2-54e748092eaa" />


The dashboard can be enabled with:

```bash
--verbose 1
```

or made quieter with:

```bash
--verbose 0
```

When verbose output is enabled, DAFT may show a rich dashboard header such as:

```text
DAFT Direction / djiNNI
```

The dashboard reports the current status of the analysis, including the number of gene trees detected, the number of completed gene trees, elapsed wall time, average time per gene tree, estimated remaining time, and total runtime.

For example, during `DAFT_Direction.py` / `djiNNI`, the dashboard may show progress for reconciling gene trees and the lineage pair currently being compared. A progress value such as `96/96` means that all 96 detected gene trees have been processed.

The line:

```text
Loading cached reconciliation for gene tree index ...
```

means that DAFT found a previously computed reconciliation result in `djiNNI_cache` and is loading it instead of recalculating that gene tree from scratch. This can make repeated runs faster when the same gene tree file and same gene tree order are used.

The field:

```text
Significant pair
```

shows the lineage pair currently being evaluated by djiNNI. For example:

```text
E; vs ((F,G),(I,C));
```

means that DAFT is comparing the direction relationship between `E;` and `((F,G),(I,C));`.

The field:

```text
Total gene trees detected
```

shows the total number of gene trees passed to the djiNNI direction step.

The field:

```text
Completed gene trees
```

shows how many gene trees have already been processed out of the total.

The field:

```text
Elapsed wall time
```

shows the real time that has passed since the current djiNNI run started.

The field:

```text
avg time/tree
```

shows the average processing time per gene tree.

The field:

```text
Remaining time
```

shows the estimated time left for the current djiNNI run.

The field:

```text
Total runtime
```

shows the estimated or final total runtime for the displayed run.

After the rich dashboard finishes, DAFT prints a final runtime summary.

```text
Total gene trees
```

is the total number of gene trees used in the direction analysis.

```text
Reconciled from scratch
```

is the number of gene trees that were newly reconciled during the current run.

```text
Loaded from cache
```

is the number of gene trees that were loaded from `djiNNI_cache` instead of being recalculated.

```text
Total elapsed wall time
```

is the final real runtime for the complete djiNNI direction analysis.

### Lineage comparison counts

After the djiNNI run, DAFT prints a section called:

```text
Lineage comparison counts
```

This section summarizes how many gene trees support each direction relationship for the significant lineage pair being tested.

For example:

```text
E; > ((F,G),(I,C));: 30
E; < ((F,G),(I,C));: 0
E; == ((F,G),(I,C));: 66
```

The `>` symbol means that the first lineage had a higher movement/support count than the second lineage for that gene tree comparison.

For example:

```text
E; > ((F,G),(I,C));: 30
```

means that 30 gene trees supported `E;` as having the higher movement/support count compared with `((F,G),(I,C));`.

The `<` symbol means that the second lineage had a higher movement/support count than the first lineage for that gene tree comparison.

For example:

```text
E; < ((F,G),(I,C));: 0
```

means that 0 gene trees supported `((F,G),(I,C));` as having the higher movement/support count compared with `E;`.

The `==` symbol means that the two lineages had equal movement/support counts for that gene tree comparison.

For example:

```text
E; == ((F,G),(I,C));: 66
```

means that 66 gene trees had equal counts between `E;` and `((F,G),(I,C));`.

In the current implementation, when the two lineages are tied under the `==` case, DAFT resolves the tie using `np.random.choice`. This means the `==` category represents gene trees where the comparison was tied or unresolved by the count comparison, and the downstream moved lineage is chosen randomly between the two lineages.

Therefore:

- `>` counts gene trees where lineage 1 had greater support.
- `<` counts gene trees where lineage 2 had greater support.
- `==` counts tied cases, which are resolved randomly using `np.random.choice`.

These counts help summarize the direction signal across gene trees before the final donor/recipient interpretation is written to the DAFT Direction output.

## 17. Troubleshooting

### `ModuleNotFoundError: No module named 'reconcILS'`

Run DAFT from the main DAFT directory:

```bash
cd DAFT
python3 DAFT_Test.py ...
```

Or pass the utility folder explicitly:

```bash
python3 DAFT_Test.py --path /path/to/DAFT_utils ...
```

By default, DAFT uses `./DAFT_utils`.

---

### `ModuleNotFoundError: No module named 'ete3'`

Install `ete3`:

```bash
pip install ete3
```

---

### `ModuleNotFoundError: No module named 'rich'`

Install `rich`:

```bash
pip install rich
```

`rich` is used for the terminal progress dashboard in the direction/transform workflow.

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

This is expected. DAFT moves outputs into folders named with the output label:

```text
DAFT_results_<output>/
DAFT_extras_<output>/
DAFT_log_<output>/
```

Check those folders after the run.

---

### Branch names are hard to interpret

Use:

```text
DAFT_results_<output>/branch_map.csv
```

to map output labels back to the original species tree.

---

### `daft-test: command not found`

Install DAFT from the repository root:

```bash
python -m pip install -e .
```

Then check the command:

```bash
which daft-test
daft-test --help
```

---

## 18. How We Recommend Reporting DAFT Results

When reporting DAFT results, include:

1. The species tree used.
2. The number of gene trees analyzed.
3. Whether `--correct 1` was used in addition to the standard run.
4. Whether sibling comparisons, avuncular comparisons, or both were used.
5. The Z-score threshold used for significance.
6. The significant lineage pairs.
7. The DAFT Direction donor/recipient inference, if direction inference was run.
8. Whether bidirectional introgression was inferred.
9. The final network, if reported.

Example language:

```text
We ran DAFT using a rooted species tree and 2,000 gene trees. 
Significance was assessed using both avuncular and sibling comparisons. 
Lineage pairs with Z <= -1.96 and sufficient supporting counts were passed to DAFT Direction/djiNNI. 
Direction inference was performed using NNI reconciliation, and the lineage with the higher movement-count was interpreted as the recipient.
```

---
