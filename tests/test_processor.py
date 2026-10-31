"""End-to-end tests on the synthetic fixture (tests/make_fixture.py), no Docker."""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import h5py
import pandas as pd
import pytest

from tests.make_fixture import MARKERS, SAMPLE_GENE, make_fixture

PARAMS = ("LABEL", "LABEL_COLUMN", "N_MARKERS", "N_TOP", "MIN_PCT", "MIN_LOGFC", "MIN_SAMPLE_FRACTION",
          "MAX_CELLS_PER_LABEL", "SEED", "THREADS")


@pytest.fixture(scope="session")
def fixture_root(tmp_path_factory):
    root = tmp_path_factory.mktemp("fixture")
    make_fixture(root, "both")
    return root


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    for k in PARAMS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv("SC_MARKERS_THREADS", "1")
    for k, v in params.items():
        monkeypatch.setenv(k, str(v))
    from processor.main import main

    return main()


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_planted_markers(fixture_root, tmp_path, monkeypatch, layout):
    assert run(fixture_root / layout, tmp_path, monkeypatch) == 0
    out = tmp_path / "cell-type"
    for name in ("markers.csv", "dotplot.png", "heatmap.png", "summary.json", "provenance.json"):
        assert (out / name).is_file(), name
    m = pd.read_csv(out / "markers.csv")
    assert list(m.columns[:7]) == ["label", "gene", "logfoldchange", "pval_adj", "pct_in", "pct_out", "rank"]
    for t, genes in MARKERS.items():
        top = m[m["label"] == t].sort_values("rank")["gene"].head(5)
        assert set(top) == set(genes), t
    assert (m.groupby("label")["rank"].min() == 1).all()
    assert (m["pval_adj"] < 0.05).all() and (m["logfoldchange"] >= 0.25).all() and (m["pct_in"] >= 0.1).all()
    # The planted sample artefact (s1 only) is dropped by the sample-aware filter.
    assert SAMPLE_GENE not in set(m["gene"])
    s = json.loads((out / "summary.json").read_text())
    assert s["n_samples"] == 4 and s["n_labels"] == 3 and s["n_dropped_inconsistent"] >= 1
    assert "Squair" in s["caveat"]
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["app"] == "sc-markers" and prov["version"] == "0.1.0"
    assert "10.1038/s41467-021-25960-2" in {c["doi"] for c in prov["citations"]}
    assert all(st["status"] == "ok" for st in prov["samples"])


def test_artefact_kept_without_sample_filter(fixture_root, tmp_path, monkeypatch):
    assert run(fixture_root / "linked", tmp_path, monkeypatch, MIN_SAMPLE_FRACTION=0) == 0
    m = pd.read_csv(tmp_path / "cell-type" / "markers.csv")
    assert SAMPLE_GENE in set(m.loc[m["label"] == "Neuron", "gene"])


def test_clusters_label_and_subsampling(fixture_root, tmp_path, monkeypatch):
    assert run(fixture_root / "linked", tmp_path, monkeypatch, LABEL_COLUMN="cluster", MAX_CELLS_PER_LABEL=100,
               N_MARKERS=3, N_TOP=2) == 0
    out = tmp_path / "cluster"
    m = pd.read_csv(out / "markers.csv")
    assert set(m["label"]) <= {"C01", "C02", "C03", "C04"}
    assert m.groupby("label").size().max() <= 3
    s = json.loads((out / "summary.json").read_text())
    assert s["n_cells_used"] == 400 and all(r["n_cells_used"] == 100 for r in s["labels"])
    # C03 and C04 split one cell type: little to tell them apart.
    assert {"C01", "C02"} <= set(m["label"])


def test_counts_in_x_are_normalized(fixture_root, tmp_path, monkeypatch):
    inp = tmp_path / "in"
    shutil.copytree(fixture_root / "linked", inp)
    for p in (inp / "cells").glob("*.h5ad"):
        with h5py.File(p, "a") as f:
            del f["X"]
            f.copy(f["layers"]["counts"], "X")
    assert run(inp, tmp_path / "out", monkeypatch) == 0
    prov = json.loads((tmp_path / "out" / "cell-type" / "provenance.json").read_text())
    assert all(s["normalized_here"] for s in prov["samples"])
    m = pd.read_csv(tmp_path / "out" / "cell-type" / "markers.csv")
    assert set(MARKERS["Neuron"]) <= set(m.loc[m["label"] == "Neuron", "gene"])


def test_missing_sample_and_bad_file(fixture_root, tmp_path, monkeypatch):
    inp = tmp_path / "in"
    shutil.copytree(fixture_root / "linked", inp)
    (inp / "cells" / "s4.h5ad").unlink()
    (inp / "cells" / "s3.h5ad").write_bytes(b"not hdf5")
    assert run(inp, tmp_path / "out", monkeypatch) == 0
    prov = json.loads((tmp_path / "out" / "cell-type" / "provenance.json").read_text())
    status = {s["sample_id"]: s["status"] for s in prov["samples"]}
    assert status == {"s1": "ok", "s2": "ok", "s3": "failed", "s4": "skipped"}


def test_symlink_escaping_input_is_skipped(fixture_root, tmp_path, monkeypatch):
    outside = tmp_path / "outside"
    outside.mkdir()
    shutil.copy(next((fixture_root / "linked" / "cells").glob("s1.h5ad")), outside / "evil.h5ad")
    inp = tmp_path / "in"
    shutil.copytree(fixture_root / "linked", inp)
    (inp / "cells" / "nested").mkdir()
    os.symlink(outside, inp / "cells" / "nested" / "link")
    assert run(inp, tmp_path / "out", monkeypatch) == 0
    prov = json.loads((tmp_path / "out" / "cell-type" / "provenance.json").read_text())
    assert any("outside INPUT_DIR" in v for v in prov["inputs"]["skipped"].values())
    assert not any("evil" in p for p in prov["inputs"]["h5ad"])


def test_bad_parameters(fixture_root, tmp_path, monkeypatch):
    assert run(fixture_root / "linked", tmp_path, monkeypatch, LABEL_COLUMN="nope") == 2
    assert run(fixture_root / "linked", tmp_path, monkeypatch, LABEL="Bad_Label") == 2
    assert run(fixture_root / "linked", tmp_path, monkeypatch, N_TOP=10, N_MARKERS=5) == 2
    assert run(fixture_root / "linked", tmp_path, monkeypatch, MIN_PCT="1.5") == 2
    assert run(fixture_root / "linked", tmp_path, monkeypatch, MAX_CELLS_PER_LABEL="x") == 2


def test_missing_inputs(tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    assert run(tmp_path / "in", tmp_path / "out", monkeypatch) == 1


def test_one_label_fails(fixture_root, tmp_path, monkeypatch):
    inp = tmp_path / "in"
    shutil.copytree(fixture_root / "linked", inp)
    p = inp / "labels" / "labels.csv"
    df = pd.read_csv(p)
    df["cell_type"] = "All"
    df.to_csv(p, index=False)
    assert run(inp, tmp_path / "out", monkeypatch) == 1


def test_handler(fixture_root, tmp_path, monkeypatch):
    for k in PARAMS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("SC_MARKERS_THREADS", "1")
    monkeypatch.setenv("N_TOP", "1")          # must be cleared by the handler
    from processor.handler import handler

    r = handler({"inputDir": str(fixture_root / "linked"), "outputDir": str(tmp_path), "LABEL": "lam"}, None)
    assert r["status"] == "ok"
    s = json.loads((tmp_path / "lam" / "summary.json").read_text())
    assert all(len(x["top"]) == 5 for x in s["labels"])
    with pytest.raises(RuntimeError):
        handler({"inputDir": str(fixture_root / "linked"), "outputDir": str(tmp_path), "LABEL_COLUMN": "x"}, None)
