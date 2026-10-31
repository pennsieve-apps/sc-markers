"""
Finding the inputs: <sample_id>.h5ad files (port `cells`, from sc-annotate or
sc-integrate) and labels.csv (port `labels`). Linked layout: INPUT_DIR/cells/,
INPUT_DIR/labels/; without links the files are found by name anywhere under
INPUT_DIR.

INPUT_DIR entries are symlinks, so the walk follows them, but only to where the
platform's own entries point: the targets of the links directly in INPUT_DIR
and in its port folders. A link further down that resolves anywhere else (a
crafted `../../etc` link inside an upstream output) is skipped.
"""
from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger("sc-markers")
SAMPLE_RE = re.compile(r"^[A-Za-z0-9_-]+$")


@dataclass
class Inputs:
    h5ad: dict[str, Path]                 # sample_id (file stem) -> path
    labels_csv: Path | None
    skipped: dict[str, str] = field(default_factory=dict)


def allowed_roots(root: Path) -> list[str]:
    """INPUT_DIR itself and where its entries (and its port folders' entries) really are."""
    roots = {os.path.realpath(root)}
    for top in root.iterdir():
        roots.add(os.path.realpath(top))
        if top.is_dir():
            for entry in top.iterdir():
                roots.add(os.path.realpath(entry))
    return sorted(roots)


def _inside(real: str, roots: list[str]) -> bool:
    return any(real == r or real.startswith(r.rstrip("/") + "/") for r in roots)


def walk(root: Path) -> tuple[list[Path], list[Path]]:
    """Every file under root (symlinks followed within the allowed roots), sorted; and the rejected paths."""
    roots = allowed_roots(root)
    found: list[Path] = []
    rejected: list[Path] = []
    seen: set[str] = set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = os.path.realpath(dirpath)
        if real in seen or not _inside(real, roots):
            if not _inside(real, roots):
                rejected.append(Path(dirpath))
            dirnames[:] = []
            continue
        seen.add(real)
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        for f in sorted(filenames):
            p = Path(dirpath) / f
            if f.startswith(".") or not p.is_file():
                continue
            (found if _inside(os.path.realpath(p), roots) else rejected).append(p)
    return found, rejected


def discover(input_dir: Path) -> Inputs:
    files, rejected = walk(input_dir)
    skipped = {str(p): "resolves outside INPUT_DIR" for p in rejected}

    def port_files(port: str) -> list[Path]:
        folder = input_dir / port
        inside = [p for p in files if folder in p.parents] if folder.is_dir() else []
        return inside or files

    h5ads = [p for p in port_files("cells") if p.suffix.lower() == ".h5ad"]
    labels = [p for p in port_files("labels") if p.name.lower() == "labels.csv"]
    by_stem: dict[str, Path] = {}
    for p in h5ads:
        if not SAMPLE_RE.match(p.stem):
            skipped[str(p)] = "file name is not a sample id ([A-Za-z0-9_-]+)"
        elif p.stem in by_stem:
            skipped[str(p)] = f"two files for sample {p.stem}; using {by_stem[p.stem]}"
        else:
            by_stem[p.stem] = p
    if len(labels) > 1:
        log.warning("%d labels.csv files; using %s", len(labels), labels[0])
    return Inputs(by_stem, labels[0] if labels else None, skipped)
