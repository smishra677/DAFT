## Running simulation experiments
(written by Sarthak Mishra, refactored by AI) 
These scripts simulate gene alignments with AliSim, infer maximum-likelihood trees with IQ-TREE, and export the inferred trees to a CSV for downstream experiments.

### Requirements

- Python 3
- IQ-TREE with AliSim support
- DendroPy: `python -m pip install dendropy`
- An input CSV containing a `gt` column of semicolon-terminated Newick trees. Quote Newick strings because they contain commas.

Run the following commands from the repository root.

### 1. Simulate alignments and infer trees

```bash
python scripts/alisim_C.py \
  --input DATA_SET/list_C.csv \
  --output DATA_SET/C_alisim_experiment \
  --iqtree ./iqtree3 \
  --model JC \
  --length 500 \
  --indel 0.004 \
  --branch-length 0.1 \
  --threads 1 \
  --seed 12345
```

This example uses the JC model, an ancestral sequence length of 500 bases, and insertion and deletion rates of 0.004 each relative to substitutions. Missing branch lengths are assigned 0.1 substitutions per site; existing lengths are retained. Indels can change the final alignment length.

Each input tree produces one simulated alignment and one inferred tree. Inference uses AliSim's true alignment directly, without realignment. The seed increases by one for each input tree.

The output directory must not already exist. It contains:

- `input_trees/`: input trees with missing branch lengths filled in.
- `alignments/`: simulated PHYLIP alignments.
- `trees/`: IQ-TREE inference outputs, including `.treefile` files.
- `logs/`: simulation and inference logs.
- `reconstructed_trees.nwk`: inferred trees, one per line.
- `manifest.csv` and `settings.json`: file mappings and run parameters.

### 2. Export inferred trees to CSV

```bash
python scripts/extract_binary_iqtree.py \
  DATA_SET/C_alisim_experiment/trees \
  --reference-csv DATA_SET/list_C.csv \
  --output DATA_SET/C_inferred.csv
```

Tree files are sorted naturally and paired with reference CSV rows by position. Use the same reference CSV and row order as the simulation input. Each target tree must have the same taxa as its reference, and each reference must have exactly two root children.

The script attempts to reroot each inferred tree using its reference root split, placing the root halfway along the matching edge. If that split is absent or rerooting fails, the current implementation resolves polytomies instead; reference-consistent rooting is therefore not guaranteed for every output tree.

The resulting CSV contains a single `gt` column. Add `--recursive` to search subdirectories or `--overwrite` to replace an existing output CSV.
