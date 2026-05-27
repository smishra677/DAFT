# Installation

DAFT requires Python 3.9–3.12.

Install from the repository root:

```bash
python -m pip install .
```

This installs the DAFT command-line tools:

```bash
daft-test
daft-direction
daft-transform
daft-excel
daft-excel-correction
daft-clean
```

Check the installation:

```bash
daft-test --help
daft-direction --help
```

## Dependencies

DAFT uses:

```text
pandas==2.2.2
numpy==1.26.4
matplotlib
igraph
ete3
openpyxl
rich
```

These are installed automatically from `setup.py`. They are also listed in `requirements.txt`.

## Example

Run DAFT Test:

```bash
daft-test \
  --sp "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));" \
  --gt "test_data.csv" \
  --output "output" \
  --sibling 1 \
  --excel 1 \
  --direction 0
```

Run DAFT Direction on selected lineage pairs:

```bash
daft-direction \
  --sp "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));" \
  --gt "test_data.csv" \
  --lineages "[(\'C;\', \'H;\')]" \
  --verbose 1 \
  --output "output"
```

## Script usage

The original scripts remain available:

```bash
python3 DAFT_Test.py --path "./DAFT_utils" ...
python3 DAFT_Direction.py --path "./DAFT_utils" ...
```

Installed commands locate `DAFT_utils` automatically.
