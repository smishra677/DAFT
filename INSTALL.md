# DAFT Installation

## Recommended installation

From the DAFT repository root:

```bash
python -m pip install -r requirements.txt
python -m pip install -e .
```

This installs DAFT with these dependencies:

```text
pandas==2.2.2
numpy==1.26.4
matplotlib
igraph
ete3
openpyxl
rich
```

The pinned `pandas` and `numpy` versions are reflected in both
`requirements.txt` and `setup.py`.

## Installed commands

After installation, these commands are available:

```bash
daft-test --help
daft-direction --help
daft-transform --help
daft-excel --help
daft-excel-correction --help
daft-clean --help
```

The installed commands automatically locate the bundled `DAFT_utils`
directory. The original script-style workflow still works:

```bash
python DAFT_Test.py ...
python DAFT_Direction.py ...
```

## Python version

This installer declares Python `>=3.9,<3.13`, because the requested
dependency set pins `numpy==1.26.4`, which is intended for Python versions
before Python 3.13.
