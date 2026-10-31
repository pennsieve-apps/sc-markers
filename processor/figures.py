"""Figures: a dot plot and a heatmap of the top markers per label (PNG, 150 dpi), drawn with scanpy."""
from __future__ import annotations

from pathlib import Path

import numpy as np

DPI = 150
HEATMAP_CELLS_PER_LABEL = 100


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "figure.facecolor": "white"})
    return plt


def dotplot(adata, genes: dict[str, list[str]], title: str, path: Path) -> Path:
    import scanpy as sc

    plt = _plt()
    dp = sc.pl.dotplot(adata, genes, groupby="label", standard_scale="var", title=title, show=False,
                       return_fig=True, colorbar_title="Mean expression\n(scaled per gene)")
    dp.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close("all")
    return path


def heatmap(adata, genes: dict[str, list[str]], title: str, path: Path, seed: int) -> Path:
    """Cells (at most HEATMAP_CELLS_PER_LABEL per label, at random) x top genes, each gene scaled 0-1."""
    import scanpy as sc

    plt = _plt()
    rng = np.random.default_rng(seed)
    codes = adata.obs["label"].cat.codes.to_numpy()
    keep = []
    for c in np.unique(codes):
        idx = np.flatnonzero(codes == c)
        keep.append(idx if len(idx) <= HEATMAP_CELLS_PER_LABEL else np.sort(rng.choice(idx, HEATMAP_CELLS_PER_LABEL,
                                                                                     replace=False)))
    sub = adata[np.concatenate(keep)].copy()
    n_genes = sum(len(v) for v in genes.values())
    axes = sc.pl.heatmap(sub, genes, groupby="label", standard_scale="var", cmap="viridis", show=False,
                         figsize=(max(6.0, 0.22 * n_genes + 3), max(4.0, 0.3 * len(genes) + 2)))
    axes["groupby_ax"].set_ylabel("")
    fig = axes["heatmap_ax"].figure
    fig.suptitle(title, fontsize=10, y=1.08)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close("all")
    return path
