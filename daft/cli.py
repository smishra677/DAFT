"""Console-script wrappers for DAFT.

These wrappers let installed users run the existing DAFT scripts from any
working directory. They do not change DAFT output formats.
"""
from __future__ import annotations

import importlib.util
import runpy
import sys
from pathlib import Path


def _has_option(flag: str) -> bool:
    """Return True when sys.argv already contains a CLI option."""
    return any(arg == flag or arg.startswith(flag + "=") for arg in sys.argv[1:])


def _daft_utils_path() -> str:
    """Locate the installed DAFT_utils directory."""
    spec = importlib.util.find_spec("DAFT_utils")
    if spec is not None:
        if spec.submodule_search_locations:
            return str(Path(next(iter(spec.submodule_search_locations))).resolve())
        if spec.origin:
            return str(Path(spec.origin).resolve().parent)

    # Source-tree fallback for editable/development use.
    candidate = Path(__file__).resolve().parents[1] / "DAFT_utils"
    return str(candidate)


def _run_module(module_name: str, needs_path: bool = False) -> None:
    """Run one of the legacy DAFT scripts as if it were called directly."""
    if needs_path and not _has_option("--path"):
        sys.argv.extend(["--path", _daft_utils_path()])
    runpy.run_module(module_name, run_name="__main__")


def daft_test() -> None:
    _run_module("DAFT_Test", needs_path=True)


def daft_direction() -> None:
    _run_module("DAFT_Direction", needs_path=True)


def daft_transform() -> None:
    _run_module("DAFT_Transform", needs_path=True)


def daft_excel() -> None:
    _run_module("DAFT_produce_excel", needs_path=False)


def daft_excel_correction() -> None:
    _run_module("DAFT_produce_excel_correction", needs_path=False)


def daft_clean() -> None:
    _run_module("clean_folder", needs_path=False)
