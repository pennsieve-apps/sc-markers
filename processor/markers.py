"""
Marker genes per label: scanpy's Wilcoxon rank-sum test (each label against all
other cells) on log-normalized expression, pooled across samples, then filtered
by effect size, detection and consistency across samples.

1. Cells. Every label's cells are taken from labels.csv; a label with more than
   MAX_CELLS_PER_LABEL cells is subsampled to that many, uniformly at random
   over its cells (so each sample keeps its share). Only the chosen rows of each
   <sample_id>.h5ad are read, a block of rows at a time; genes are those present
   in every sample.
2. Test. scanpy.tl.rank_genes_groups(method="wilcoxon", tie_correct=True,
   pts=True), Benjamini-Hochberg per label, on the genes detected in at least
   MIN_PCT of the cells of some label (the others cannot pass the MIN_PCT filter).
3. Sample-aware filter. Pooling cells from several donors makes the test treat
   cells, not samples, as replicates: a gene that is high in one label in a
   single donor can reach a tiny p-value. So, within each sample that has at
   least MIN_CELLS_SAMPLE cells of the label and of the rest, the label's log2
   fold change over the rest of that sample's cells is computed; a marker must
   reach MIN_LOGFC in at least MIN_SAMPLE_FRACTION of those samples.

A marker: p_adj < 0.05, log2 fold change >= MIN_LOGFC, detected in >= MIN_PCT of
the label's cells, and the sample-aware filter. Ranked by the Wilcoxon z-score.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

log = logging.getLogger("sc-markers")
CHUNK_ROWS = 20000
P_ADJ = 0.05
MIN_CELLS_LABEL = 10          # a label with fewer cells is not tested
MIN_CELLS_SAMPLE = 10         # a sample takes part in the consistency check for a label with this many of its cells (and of the rest)
TARGET_SUM = 1e4


@dataclass
class SampleBlock:
    sample_id: str
    genes: list[str]
    X: sp.csr_matrix              # chosen rows x genes, log-normalized float32
    obs_ids: list[str]
    n_obs_file: int
    normalized_here: bool


def _read_index(group) -> np.ndarray:
    from anndata.io import read_elem

    key = group.attrs.get("_index", "_index")
    return np.asarray(read_elem(group[key])).astype(str)


def _iter_rows(elem, n_obs: int, n_vars: int, chunk: int):
    """Yield (start, stop, csr block) over the rows of a CSR group or a dense dataset."""
    import h5py

    enc = elem.attrs.get("encoding-type", "")
    if isinstance(elem, h5py.Group) and enc == "csr_matrix":
        indptr = elem["indptr"][:]
        for start in range(0, n_obs, chunk):
            stop = min(start + chunk, n_obs)
            a, b = int(indptr[start]), int(indptr[stop])
            yield start, stop, sp.csr_matrix((elem["data"][a:b], elem["indices"][a:b], indptr[start:stop + 1] - a),
                                             shape=(stop - start, n_vars))
    elif isinstance(elem, h5py.Dataset):
        for start in range(0, n_obs, chunk):
            stop = min(start + chunk, n_obs)
            yield start, stop, sp.csr_matrix(elem[start:stop])
    else:
        raise ValueError(f"X is stored as {enc or type(elem).__name__}; expected CSR or dense")


def read_sample(path: Path, sample_id: str, wanted: pd.Index, chunk: int = CHUNK_ROWS) -> SampleBlock:
    """The rows of X whose obs_names are in `wanted`, log-normalized (X as written by sc-qc and later steps;
    integer counts, as sc-prep writes them, are normalized to TARGET_SUM and log1p'd here)."""
    import h5py

    with h5py.File(path, "r") as f:
        obs_names = _read_index(f["obs"])
        genes = _read_index(f["var"]).tolist()
        take = np.flatnonzero(pd.Index(obs_names).isin(wanted))
        blocks, is_int = [], True
        for start, stop, block in _iter_rows(f["X"], len(obs_names), len(genes), chunk):
            rows = take[(take >= start) & (take < stop)] - start
            if rows.size == 0:
                continue
            b = block[rows].astype(np.float32)
            if b.data.size and b.data.min() < 0:
                raise ValueError("X has negative values (scaled data?); expected log-normalized expression")
            is_int = is_int and bool(np.all(np.mod(b.data, 1) == 0))
            blocks.append(b)
    X = sp.vstack(blocks, format="csr") if blocks else sp.csr_matrix((0, len(genes)), dtype=np.float32)
    normalized_here = bool(X.nnz) and is_int
    if normalized_here:
        totals = np.asarray(X.sum(axis=1)).ravel()
        X = sp.csr_matrix(sp.diags(TARGET_SUM / np.maximum(totals, 1)) @ X, dtype=np.float32)
        X.data = np.log1p(X.data)
    return SampleBlock(sample_id, genes, X, obs_names[take].tolist(), len(obs_names), normalized_here)


def choose_cells(labels: pd.DataFrame, max_per_label: int, seed: int) -> pd.DataFrame:
    """At most max_per_label cells of each label, uniformly at random over the label's cells."""
    rng = np.random.default_rng(seed)
    keep = []
    for _, g in labels.groupby("label", sort=True, observed=True):
        idx = g.index.to_numpy()
        if len(idx) > max_per_label:
            idx = np.sort(rng.choice(idx, size=max_per_label, replace=False))
        keep.append(idx)
    return labels.loc[np.concatenate(keep)] if keep else labels.iloc[0:0]


def sample_consistency(X: sp.csr_matrix, sample: np.ndarray, label: np.ndarray, labels: list[str],
                       min_cells: int, min_logfc: float) -> tuple[np.ndarray, np.ndarray]:
    """Per label x gene: in how many samples the label's log2 fold change over the rest of that sample's cells
    (computed as scanpy does, from the means of expm1(X)) is at least min_logfc (n_up); and per label how many
    samples had >= min_cells of the label and of the rest (n_tested)."""
    samples = sorted(set(sample))
    s_code = pd.Categorical(sample, categories=samples).codes
    l_code = pd.Categorical(label, categories=labels).codes
    n_s, n_l = len(samples), len(labels)
    grp = s_code * n_l + l_code
    ind = sp.csr_matrix((np.ones(len(grp)), (grp, np.arange(len(grp)))), shape=(n_s * n_l, len(grp)))
    E = X.copy()
    E.data = np.expm1(E.data)
    sums = np.asarray((ind @ E).todense()).reshape(n_s, n_l, -1)            # sample x label x gene
    counts = np.bincount(grp, minlength=n_s * n_l).reshape(n_s, n_l).astype(np.float64)
    tot_sum, tot_n = sums.sum(axis=1), counts.sum(axis=1)                    # sample x gene, sample
    n_up = np.zeros((n_l, X.shape[1]), dtype=np.int32)
    n_tested = np.zeros(n_l, dtype=np.int32)
    for li in range(n_l):
        for si in range(n_s):
            n_in, n_out = counts[si, li], tot_n[si] - counts[si, li]
            if n_in < min_cells or n_out < min_cells:
                continue
            n_tested[li] += 1
            lfc = np.log2((sums[si, li] / n_in + 1e-9) / ((tot_sum[si] - sums[si, li]) / n_out + 1e-9))
            n_up[li] += lfc >= min_logfc
    return n_up, n_tested


def find_markers(adata, params: dict, n_threads: int) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """adata: cells x genes, log-normalized X, obs 'label' (categorical) and 'sample_id'.
    Returns (markers, per-label table, number of genes tested)."""
    import scanpy as sc

    sc.settings.n_jobs = n_threads
    labels = list(adata.obs["label"].cat.categories)
    lab_code = adata.obs["label"].cat.codes.to_numpy()
    # Genes detected in >= MIN_PCT of some label's cells; the others can never pass the detection filter.
    ind = sp.csr_matrix((np.ones(adata.n_obs), (lab_code, np.arange(adata.n_obs))), shape=(len(labels), adata.n_obs))
    det = np.asarray((ind @ (adata.X > 0).astype(np.float32)).todense())
    frac = det / np.bincount(lab_code, minlength=len(labels))[:, None]
    gene_mask = (frac >= params["MIN_PCT"]).any(axis=0)
    ad = adata[:, gene_mask].copy()
    log.info("testing %d of %d genes (detected in >= %.0f%% of some label's cells)", ad.n_vars, adata.n_vars,
             100 * params["MIN_PCT"])
    sc.tl.rank_genes_groups(ad, "label", method="wilcoxon", tie_correct=True, pts=True, use_raw=False)
    res = sc.get.rank_genes_groups_df(ad, group=None)
    res = res.rename(columns={"group": "label", "names": "gene", "scores": "score", "logfoldchanges": "logfoldchange",
                              "pvals": "pval", "pvals_adj": "pval_adj", "pct_nz_group": "pct_in",
                              "pct_nz_reference": "pct_out"})
    res["label"] = res["label"].astype(str)

    n_up, n_tested = sample_consistency(ad.X, ad.obs["sample_id"].astype(str).to_numpy(),
                                        ad.obs["label"].astype(str).to_numpy(), labels, MIN_CELLS_SAMPLE,
                                        params["MIN_LOGFC"])
    gpos = pd.Index(ad.var_names).get_indexer(res["gene"])
    lpos = pd.Index(labels).get_indexer(res["label"])
    res["n_samples_up"] = n_up[lpos, gpos]
    res["n_samples"] = n_tested[lpos]
    need = np.ceil(params["MIN_SAMPLE_FRACTION"] * res["n_samples"].to_numpy() - 1e-9)
    res["consistent"] = res["n_samples_up"].to_numpy() >= need
    passed = ((res["pval_adj"] < P_ADJ) & (res["logfoldchange"] >= params["MIN_LOGFC"])
              & (res["pct_in"] >= params["MIN_PCT"]))
    res["passed_pooled"] = passed
    markers = res[passed & res["consistent"]].sort_values(["label", "score"], ascending=[True, False])
    markers = markers.groupby("label", sort=False).head(params["N_MARKERS"]).copy()
    markers["rank"] = markers.groupby("label", sort=False).cumcount() + 1
    per_label = pd.DataFrame({
        "label": labels,
        "n_pooled": [int(passed[res["label"] == lab].sum()) for lab in labels],
        "n_dropped_inconsistent": [int((passed & ~res["consistent"])[res["label"] == lab].sum()) for lab in labels],
        "n_samples_checked": n_tested,
    })
    cols = ["label", "gene", "logfoldchange", "pval_adj", "pct_in", "pct_out", "rank", "score", "n_samples_up",
            "n_samples"]
    return markers[cols].reset_index(drop=True), per_label, int(ad.n_vars)


def top_genes(markers: pd.DataFrame, n_top: int) -> dict[str, list[str]]:
    """Top n_top markers per label for the figures; a gene is shown once, under the first label that ranks it."""
    out: dict[str, list[str]] = {}
    used: set[str] = set()
    for lab, g in markers.groupby("label", sort=False):
        genes = [x for x in g.sort_values("rank")["gene"] if x not in used][:n_top]
        used.update(genes)
        if genes:
            out[str(lab)] = genes
    return out
