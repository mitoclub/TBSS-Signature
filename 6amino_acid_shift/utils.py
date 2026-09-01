"""Shared utilities for MIDORI2 amino-acid ratio analysis."""

from __future__ import annotations

import hashlib
from pathlib import Path


TARGET_GENES = (
    "CO1", "CO2", "A8", "A6", "CO3",
    "ND3", "ND4L", "ND4", "ND5", "Cytb",
)


def sha256(path: Path) -> str:
    """Calculate SHA256 checksum of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = ["TARGET_GENES", "sha256"]
