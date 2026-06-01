"""DAFT utility package.

This file makes DAFT_utils importable as a package while preserving DAFT's
legacy script-style imports.
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_RECONCILS = _HERE / "reconcILS"

# Preserve legacy imports used inside DAFT_essential.py:
#   from utils_reconcILS import *
#   from reconcILS import *
for _path in (_HERE, _RECONCILS):
    _path_str = str(_path)
    if _path_str not in sys.path:
        sys.path.insert(0, _path_str)


def __getattr__(name: str):
    """Lazy export so importing DAFT_utils does not immediately load everything."""
    if name == "daft_essential":
        from .DAFT_essential import daft_essential
        return daft_essential

    elif name == "daft_validate":
        from .DAFT_validate import daft_validate
        return daft_validate
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["daft_essential","daft_validate"]