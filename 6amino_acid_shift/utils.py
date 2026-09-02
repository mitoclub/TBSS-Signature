"""Shared utilities for MIDORI2 amino-acid ratio analysis."""

from __future__ import annotations

import hashlib
from pathlib import Path
import sys

_REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
if str(_REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPOSITORY_ROOT))

from mtdna import MAJOR_ARC_HEAVY_STRAND

# COX1 -> CytB in canonical rCRS order, defined once in ../mtdna.py.
TARGET_GENES = MAJOR_ARC_HEAVY_STRAND


def sha256(path: Path) -> str:
    """Calculate SHA256 checksum of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = ["TARGET_GENES", "sha256"]
