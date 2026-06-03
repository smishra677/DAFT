# Installing DAFT

This page gives the recommended installation steps for DAFT. The main README gives the full usage guide and output-format documentation.

---

## Requirements

DAFT requires **Python 3.9–3.12**.

The Python dependencies are listed in both `setup.py` and `requirements.txt`:

| Package | Version / requirement | Used for |
|---|---:|---|
| `pandas` | `2.2.2` | Reading gene-tree tables and writing tabular outputs |
| `numpy` | `1.26.4` | Numerical calculations and reproducible random tie-breaking |
| `matplotlib` | latest compatible version | Plotting / visualization support |
| `igraph` | latest compatible version | Network construction and graph handling |
| `ete3` | latest compatible version | Parsing, rooting, and traversing Newick trees |
| `openpyxl` | latest compatible version | Excel output support |
| `rich` | latest compatible version | Terminal progress display |

---

## Recommended installation

From the DAFT repository root, install DAFT with `pip`:

```bash
python -m pip install .
```

This installs DAFT and its Python dependencies. After installation, the main commands are available from the command line.

| Command | Purpose |
|---|---|
| `daft-test` | Run DAFT Test |
| `daft-direction` | Run DAFT Direction |
| `daft-transform` | Run the DAFT transform / djiNNI workflow |
| `daft-excel` | Produce Excel output from DAFT results |
| `daft-excel-correction` | Produce Excel output for corrected DAFT results |
| `daft-clean` | Clean DAFT output folders |

Check that the installation worked:

```bash
daft-test --help
daft-direction --help
daft-transform --help
```

The installed commands automatically locate `DAFT_utils`, so users normally do **not** need to pass `--path` when using the installed command-line tools.

---

## Optional: use a virtual environment

A virtual environment is recommended if you do not want DAFT dependencies to affect other Python projects.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

Then install DAFT:

```powershell
python -m pip install .
```

---

## Optional: editable/development installation

If you are editing the DAFT source code and want command-line changes to take effect without reinstalling each time, use editable mode:

```bash
python -m pip install -e .
```

This is the recommended installation mode for development.

---

## Optional: install dependencies separately

Usually this is not necessary, because `python -m pip install .` installs the dependencies from `setup.py`.

If you want to install dependencies first, run:

```bash
python -m pip install -r requirements.txt
python -m pip install .
```

---

## Example test run

After installing DAFT, you can test the installation using the included `test_data.csv` file.

```bash
daft-test \
  --sp "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));" \
  --gt "test_data.csv" \
  --output "output" \
  --sibling 1 \
  --excel 1 \
  --direction 0
```

To run DAFT Direction on selected lineage pairs:

```bash
daft-direction \
  --sp "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));" \
  --gt "test_data.csv" \
  --lineages "[('C;', 'H;')]" \
  --verbose 1 \
  --output "output"
```

---

## Running the original scripts directly

The installed command-line tools are recommended. However, the original scripts can still be run directly from the repository root:

```bash
python DAFT_Test.py --path "./DAFT_utils" ...
python DAFT_Direction.py --path "./DAFT_utils" ...
python DAFT_Transform.py --path "./DAFT_utils" ...
```

When running the original scripts from outside the repository root, pass the absolute path to `DAFT_utils`:

```bash
python /path/to/DAFT/DAFT_Test.py --path /path/to/DAFT/DAFT_utils ...
```

---

## Troubleshooting

### `daft-test: command not found`

DAFT may not be installed in the active Python environment. From the repository root, run:

```bash
python -m pip install .
```

Then check:

```bash
daft-test --help
```

If you are using a virtual environment, make sure it is activated before running DAFT.

---

### `ModuleNotFoundError: No module named 'DAFT_utils'`

This usually means you are running the original Python scripts directly and DAFT cannot find the utility folder.

Preferred fix: use the installed command-line tools:

```bash
daft-test --help
```

Alternative fix: pass `--path` explicitly:

```bash
python DAFT_Test.py --path /path/to/DAFT/DAFT_utils ...
```

---

### `ModuleNotFoundError` for `ete3`, `igraph`, `rich`, `openpyxl`, `pandas`, or `numpy`

Install DAFT from the repository root so dependencies are installed automatically:

```bash
python -m pip install .
```

Or install the requirements directly:

```bash
python -m pip install -r requirements.txt
```

If installation still fails, upgrade `pip` first:

```bash
python -m pip install --upgrade pip
```

---

### `Permission denied` during installation

Use a virtual environment instead of installing into the system Python:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
```

---

## Updating DAFT

After pulling or downloading a newer version of DAFT, reinstall from the repository root:

```bash
python -m pip install --upgrade .
```

For an editable/development installation, run:

```bash
python -m pip install -e .
```

---

## Uninstalling DAFT

```bash
python -m pip uninstall daft-introgression
```
