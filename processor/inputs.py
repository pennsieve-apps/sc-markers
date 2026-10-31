"""
Finding the inputs: <sample_id>.h5ad files (port `cells`, from sc-annotate or sc-integrate) and
labels.csv (port `labels`). Linked layout: INPUT_DIR/cells/, INPUT_DIR/labels/; without links the
files are found by name anywhere under INPUT_DIR. The files come from safe_inputs.InputDir.walk(),
which applies the platform's symlink rule (conventions.md).
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

from processor.safe_inputs import InputDir

log = logging.getLogger("sc-markers")
SAMPLE_RE = re.compile(r"^[A-Za-z0-9_-]+$")


@dataclass
class Inputs:
    dir: InputDir
    h5ad: dict[str, Path]                 # sample_id (file stem) -> path
    labels_csv: Path | None
    skipped: dict[str, str] = field(default_factory=dict)


def discover(input_dir: Path) -> Inputs:
    d = InputDir(input_dir)
    files = d.walk()

    def port_files(port: str) -> list[Path]:
        folder = input_dir / port
        inside = [p for p in files if folder in p.parents]
        return inside or files

    h5ads = [p for p in port_files("cells") if p.suffix.lower() == ".h5ad"]
    labels = [p for p in port_files("labels") if p.name.lower() == "labels.csv"]
    by_stem: dict[str, Path] = {}
    skipped: dict[str, str] = {}
    for p in h5ads:
        if not SAMPLE_RE.match(p.stem):
            skipped[d.rel(p)] = "file name is not a sample id ([A-Za-z0-9_-]+)"
        elif p.stem in by_stem:
            skipped[d.rel(p)] = f"two files for sample {p.stem}; using {d.rel(by_stem[p.stem])}"
        else:
            by_stem[p.stem] = p
    if len(labels) > 1:
        log.warning("%d labels.csv files; using %s", len(labels), d.rel(labels[0]))
    return Inputs(d, by_stem, labels[0] if labels else None, skipped)
