"""
A synthetic fixture that follows the single-cell conventions (sc-annotate's `cells` and `labels`
outputs), with planted markers:

  4 samples (s1, s2 in condition `control`, s3, s4 in `pain`), 400 nuclei each, 300 genes.
  3 cell types: Neuron (15 %), Satellite_glia (35 %), Fibroblast (50 %), each with five marker
  genes 4x higher than elsewhere; clusters C01 (Neuron), C02 (Satellite_glia), C03 and C04
  (Fibroblast, split).
  Planted sample artefact: SAMPLE_GENE is 12x higher in s1's Neuron nuclei only, so a pooled
  test can call it a Neuron marker; the sample-aware filter must drop it (1 of 4 samples).

Writes INPUT_DIR layouts: linked (cells/, labels/), flat (everything side by side), and platform:
the real link layout, platform/cells/ and platform/labels/ real folders whose entries are symlinks
into a sibling upstream output folder (upstream/cells/, upstream/labels/).

    python -m tests.make_fixture data/fixture --layout both
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

SAMPLES = [("s1", "control"), ("s2", "control"), ("s3", "pain"), ("s4", "pain")]
TYPES = ["Neuron", "Satellite_glia", "Fibroblast"]
SHARES = [0.15, 0.35, 0.50]
N_GENES = 300
N_CELLS = 400
MARKERS = {"Neuron": ["SNAP25", "NEFL", "PRPH", "SCN10A", "CALCA"],
           "Satellite_glia": ["FABP7", "APOE", "SOX10", "S100B", "PLP1"],
           "Fibroblast": ["DCN", "LUM", "COL1A1", "PDGFRA", "APOD"]}
SAMPLE_GENE = "SAMPLE_GENE"
MT_GENES = ["MT-CO1", "MT-ND1", "MT-ATP6"]


def gene_names() -> list[str]:
    named = [g for gs in MARKERS.values() for g in gs] + [SAMPLE_GENE] + MT_GENES
    return named + [f"G{i:03d}" for i in range(N_GENES - len(named))]


BASE = np.exp(np.random.default_rng(0).normal(np.log(0.3), 1.0, N_GENES))


def make_sample(sample_id: str, condition: str, rng: np.random.Generator):
    import anndata as ad

    genes = gene_names()
    types = rng.choice(TYPES, size=N_CELLS, p=SHARES)
    counts = np.zeros((N_CELLS, N_GENES), dtype=np.int32)
    sample_factor = np.exp(rng.normal(0.0, 0.15, N_GENES))
    for t in TYPES:
        rows = np.flatnonzero(types == t)
        mu = BASE * sample_factor
        mu[[genes.index(g) for g in MARKERS[t]]] *= 4.0
        if t == "Neuron" and sample_id == "s1":
            mu[genes.index(SAMPLE_GENE)] *= 12.0
        depth = rng.gamma(4.0, 0.25, size=(len(rows), 1))
        counts[rows] = rng.poisson(mu[None, :] * depth * 6)
    cluster = np.where(types == "Neuron", "C01", np.where(types == "Satellite_glia", "C02",
                                                           np.where(rng.random(N_CELLS) < 0.5, "C03", "C04")))
    barcodes = [f"AAAC{i:08d}-1" for i in range(N_CELLS)]
    obs_names = [f"{sample_id}:{b}" for b in barcodes]
    X = sp.csr_matrix(counts, dtype=np.float32)
    totals = np.asarray(X.sum(axis=1)).ravel()
    norm = sp.csr_matrix(sp.diags(1e4 / np.maximum(totals, 1)) @ X, dtype=np.float32)
    norm.data = np.log1p(norm.data)
    mt = np.array([g in MT_GENES for g in genes])
    obs = pd.DataFrame({
        "sample_id": pd.Categorical([sample_id] * N_CELLS),
        "total_counts": totals.astype(np.float32),
        "n_genes_by_counts": np.asarray((X > 0).sum(axis=1)).ravel().astype(np.int32),
        "pct_counts_mt": (100 * np.asarray(X[:, mt].sum(axis=1)).ravel() / np.maximum(totals, 1)).astype(np.float32),
        "pct_counts_ribo": np.float32(0.0),
        "doublet_score": rng.random(N_CELLS).astype(np.float32) * 0.2,
        "predicted_doublet": False,
        "cluster": pd.Categorical(cluster),
        "cell_type": pd.Categorical(types),
        "cell_type_score": rng.uniform(0.6, 1.0, N_CELLS).astype(np.float32),
        "cell_type_source": pd.Categorical(["markers"] * N_CELLS),
        "condition": pd.Categorical([condition] * N_CELLS),
        "donor": pd.Categorical([f"d{sample_id[1:]}"] * N_CELLS),
    }, index=obs_names)
    var = pd.DataFrame({"gene_ids": [f"ENSG{i:011d}" for i in range(N_GENES)], "feature_types": "Gene Expression",
                        "mt": mt, "ribo": False}, index=genes)
    adata = ad.AnnData(X=norm, obs=obs, var=var)
    adata.layers["counts"] = sp.csr_matrix(counts, dtype=np.float32)
    emb = rng.normal(size=(N_CELLS, 10)).astype(np.float32)
    adata.obsm["X_pca"] = emb
    adata.obsm["X_integrated"] = emb
    adata.obsm["X_umap"] = rng.normal(size=(N_CELLS, 2)).astype(np.float32)
    adata.uns["pennsieve"] = {"source": f"fixture/{sample_id}", "chemistry": "flex", "organism": "human",
                              "prep_version": "fixture"}
    labels = pd.DataFrame({"sample_id": sample_id, "obs_id": obs_names, "cluster": cluster, "cell_type": types,
                           "cell_type_score": obs["cell_type_score"].to_numpy()})
    return adata, labels


def make_fixture(root: Path, layout: str = "both", seed: int = 7) -> dict[str, Path]:
    root = Path(root)
    rng = np.random.default_rng(seed)
    made = [(sid, *make_sample(sid, cond, rng)) for sid, cond in SAMPLES]
    labels = pd.concat([m[2] for m in made], ignore_index=True)
    out: dict[str, Path] = {}
    for name in (["linked", "flat", "platform"] if layout == "both" else [layout]):
        base = root / ("upstream" if name == "platform" else name)
        if base.exists():
            shutil.rmtree(base)
        cells = base / "cells" if name != "flat" else base
        lab = base / "labels" if name != "flat" else base
        for d in (cells, lab):
            d.mkdir(parents=True, exist_ok=True)
        for sid, adata, _ in made:
            adata.write_h5ad(cells / f"{sid}.h5ad")
        labels.to_csv(lab / "labels.csv", index=False)
        out[name] = base
        if name == "platform":
            link_root = root / "platform"
            if link_root.exists():
                shutil.rmtree(link_root)
            for port in ("cells", "labels"):
                (link_root / port).mkdir(parents=True)
                for f in sorted((base / port).iterdir()):
                    (link_root / port / f.name).symlink_to(Path("..", "..", "upstream", port, f.name))
            out[name] = link_root
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out", type=Path)
    parser.add_argument("--layout", choices=["linked", "flat", "platform", "both"], default="both")
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args(argv)
    for name, path in make_fixture(args.out, args.layout, args.seed).items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
