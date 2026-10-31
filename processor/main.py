"""
sc-markers: which genes mark each cell type (or cluster)?

Reads <sample_id>.h5ad files (port `cells`, log-normalized X, from sc-annotate
or sc-integrate) and labels.csv (port `labels`), and writes, all under
OUTPUT_DIR/<LABEL>/ so several runs (cell types, clusters) can be linked into
one report folder:

  markers.csv      label, gene, logfoldchange, pval_adj, pct_in, pct_out, rank, score, n_samples_up, n_samples
  dotplot.png      top N_TOP markers per label: mean expression (scaled per gene) and fraction of cells
  heatmap.png      the same genes in single cells (at most 100 per label)
  summary.json, provenance.json

See processor/markers.py for the method.
"""
from __future__ import annotations

import json
import logging
import math
import os
import platform
import re
import sys
import time
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("sc-markers")
log.setLevel(logging.INFO)   # the Lambda runtime sets the root logger to WARNING first
for noisy in ("matplotlib", "PIL", "fontTools", "h5py", "numba", "scanpy", "anndata"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "sc-markers"
REPO_ROOT = Path(__file__).resolve().parent.parent
LABEL_RE = re.compile(r"^[a-z0-9-]+$")
MAX_LABELS_CSV_BYTES = 256 * 1024 ** 2

# The methods this step uses, in the conventions' short form (verified with scripts/verify_citations.py;
# sc-report renders the full reference by DOI).
CITATIONS = [
    {"doi": "10.1186/s13059-017-1382-0", "pmid": "29409532", "text": "Wolf FA et al. Genome Biol 2018"},
    {"doi": "10.21105/joss.04371", "text": "Virshup I et al. J Open Source Softw 2024"},
    {"doi": "10.2307/3001968", "text": "Wilcoxon F. Biometrics Bull 1945"},
    {"doi": "10.1214/aoms/1177730491", "text": "Mann HB et al. Ann Math Stat 1947"},
    {"doi": "10.1111/j.2517-6161.1995.tb02031.x", "text": "Benjamini Y et al. J R Stat Soc Series B 1995"},
    {"doi": "10.1186/s13059-024-03183-0", "pmid": "38409056", "text": "Pullin JM et al. Genome Biol 2024"},
    {"doi": "10.1038/s41467-021-25960-2", "pmid": "34584091", "text": "Squair JW et al. Nat Commun 2021"},
    {"doi": "10.1038/s41467-021-21038-1", "pmid": "33531494", "text": "Zimmerman KD et al. Nat Commun 2021"},
    {"doi": "10.1016/j.cels.2019.07.012", "pmid": "31521605", "text": "Zhang JM et al. Cell Syst 2019"},
    {"doi": "10.1038/s41576-023-00586-w", "pmid": "37002403", "text": "Heumos L et al. Nat Rev Genet 2023"},
    {"doi": "10.1109/MCSE.2007.55", "text": "Hunter JD. Comput Sci Eng 2007"},
]

DEFAULTS = {
    "LABEL": "",
    "LABEL_COLUMN": "cell_type",
    "N_MARKERS": "50",
    "N_TOP": "5",
    "MIN_PCT": "0.1",
    "MIN_LOGFC": "0.25",
    "MIN_SAMPLE_FRACTION": "0.5",
    "MAX_CELLS_PER_LABEL": "2000",
    "SEED": "42",
    "THREADS": "0",
}
PARAMETERS = tuple(DEFAULTS)


class ParameterError(ValueError):
    pass


def read_parameters(env=os.environ) -> dict:
    raw = {k: (str(env.get(k)).strip() if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    p: dict = {"LABEL_COLUMN": raw["LABEL_COLUMN"]}
    if not p["LABEL_COLUMN"] or p["LABEL_COLUMN"] in ("sample_id", "obs_id"):
        raise ParameterError("LABEL_COLUMN must name a label column of labels.csv (e.g. cell_type, cluster)")
    label = raw["LABEL"] or re.sub(r"[^a-z0-9-]+", "-", p["LABEL_COLUMN"].lower().replace("_", "-")).strip("-")
    if not LABEL_RE.match(label):
        raise ParameterError(f"LABEL must match [a-z0-9-]+, got {label!r}")
    p["LABEL"] = label
    try:
        for k in ("N_MARKERS", "N_TOP", "MAX_CELLS_PER_LABEL", "SEED", "THREADS"):
            p[k] = int(raw[k])
        for k in ("MIN_PCT", "MIN_LOGFC", "MIN_SAMPLE_FRACTION"):
            p[k] = float(raw[k])
    except ValueError as exc:
        raise ParameterError(f"bad numeric parameter: {exc}") from exc
    if p["N_MARKERS"] < 1 or p["N_TOP"] < 1 or p["N_TOP"] > p["N_MARKERS"]:
        raise ParameterError("need 1 <= N_TOP <= N_MARKERS")
    if p["MAX_CELLS_PER_LABEL"] < 50:
        raise ParameterError("MAX_CELLS_PER_LABEL must be >= 50")
    if not 0 <= p["MIN_PCT"] < 1 or not 0 <= p["MIN_SAMPLE_FRACTION"] <= 1:
        raise ParameterError("MIN_PCT must be in [0, 1) and MIN_SAMPLE_FRACTION in [0, 1]")
    if p["MIN_LOGFC"] < 0 or p["THREADS"] < 0:
        raise ParameterError("MIN_LOGFC and THREADS must be >= 0")
    return p


def _clean_json(obj):
    import numpy as np

    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        obj = float(obj)
        return None if math.isnan(obj) or math.isinf(obj) else obj
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, dict):
        return {str(k): _clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean_json(v) for v in obj]
    return obj


def write_json(path: Path, data) -> None:
    with open(path, "w") as fh:
        json.dump(_clean_json(data), fh, indent=2)
        fh.write("\n")


def tool_versions() -> dict:
    out = {"python": platform.python_version()}
    for pkg in ("scanpy", "anndata", "h5py", "numpy", "scipy", "pandas", "matplotlib"):
        try:
            out[pkg] = version(pkg)
        except PackageNotFoundError:
            out[pkg] = None
    return out


def peak_rss_mb():
    import resource

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(peak / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)


def cgroup_peak_mb():
    try:
        return round(int(open("/sys/fs/cgroup/memory.peak").read()) / 1024 / 1024, 1)
    except (OSError, ValueError):
        return None


def app_version() -> str:
    try:
        return (REPO_ROOT / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def main() -> int:
    started = datetime.now(timezone.utc)
    t0 = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba")
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("%s", exc)
        return 2
    n_threads = int(os.environ.get("SC_MARKERS_THREADS") or params["THREADS"] or 1)

    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    from processor import figures as F
    from processor import markers as M
    from processor.inputs import discover

    from processor.safe_inputs import UnsafeInput

    found = discover(input_dir)
    for path, why in found.skipped.items():
        log.warning("skipped %s: %s", path, why)
    if found.labels_csv is None:
        log.error("no labels.csv (port labels) under INPUT_DIR")
        return 1
    if not found.h5ad:
        log.error("no <sample_id>.h5ad files (port cells) under INPUT_DIR")
        return 1
    lc = params["LABEL_COLUMN"]
    try:
        header = found.dir.read_csv(found.labels_csv, MAX_LABELS_CSV_BYTES, nrows=0).columns
        missing = [c for c in ("sample_id", "obs_id") if c not in header]
        if missing:
            log.error("labels.csv lacks %s", ", ".join(missing))
            return 1
        if lc not in header:
            log.error("LABEL_COLUMN %r is not a column of labels.csv (columns: %s)", lc, ", ".join(header))
            return 2
        labels = found.dir.read_csv(found.labels_csv, MAX_LABELS_CSV_BYTES, usecols=["sample_id", "obs_id", lc],
                                    dtype=str, keep_default_na=False)
    except (UnsafeInput, ValueError) as exc:
        log.error("cannot read labels.csv: %s", exc)
        return 1
    dup = int(labels.duplicated(["sample_id", "obs_id"]).sum())
    if dup:
        log.error("labels.csv has %d duplicated (sample_id, obs_id) rows; each cell must be listed once", dup)
        return 1
    labels = labels.rename(columns={lc: "label"})
    labels["label"] = labels["label"].str.strip()
    n_unlabelled = int((labels["label"].isin(["", "NA", "nan", "NaN"])).sum())
    labels = labels[~labels["label"].isin(["", "NA", "nan", "NaN"])].reset_index(drop=True)
    notes: list[str] = []
    if n_unlabelled:
        notes.append(f"{n_unlabelled} cells without a {lc} were left out.")

    # ---- which cells -----------------------------------------------------------------------------
    statuses: dict[str, dict] = {}
    for sid in sorted(set(labels["sample_id"]) - set(found.h5ad)):
        statuses[sid] = {"sample_id": sid, "status": "skipped", "reason": "in labels.csv but no <sample_id>.h5ad"}
        log.warning("%s: in labels.csv but no %s.h5ad; its cells are left out", sid, sid)
    labels = labels[labels["sample_id"].isin(found.h5ad)].reset_index(drop=True)
    if labels.empty:
        log.error("no labelled cell belongs to a sample with an .h5ad (labels.csv sample_id = file stem?)")
        return 1
    n_per_label = labels["label"].value_counts()
    chosen = M.choose_cells(labels, params["MAX_CELLS_PER_LABEL"], params["SEED"])
    log.info("%d labelled cells in %d samples, %d labels; %d cells chosen (MAX_CELLS_PER_LABEL=%d)", len(labels),
             labels["sample_id"].nunique(), len(n_per_label), len(chosen), params["MAX_CELLS_PER_LABEL"])

    # ---- read ------------------------------------------------------------------------------------
    t_read = time.monotonic()
    blocks: list = []
    for sid, g in chosen.groupby("sample_id", sort=True):
        t = time.monotonic()
        try:
            b = M.read_sample(found.dir, found.h5ad[sid], sid, pd.Index(g["obs_id"]))
        except Exception as exc:   # one bad file is skipped, not fatal
            statuses[sid] = {"sample_id": sid, "status": "failed", "error": f"{type(exc).__name__}: {exc}"}
            log.error("%s: failed: %s", sid, exc)
            continue
        statuses[sid] = {"sample_id": sid, "status": "ok", "h5ad": found.dir.rel(found.h5ad[sid]),
                         "n_obs_file": b.n_obs_file, "n_chosen": len(g), "n_read": b.X.shape[0],
                         "n_genes": len(b.genes), "normalized_here": b.normalized_here,
                         "seconds": round(time.monotonic() - t, 2)}
        if b.X.shape[0] < len(g):
            log.warning("%s: %d of %d chosen cells not in the h5ad (obs_id = obs_names?)", sid, len(g) - b.X.shape[0],
                        len(g))
        if b.normalized_here:
            log.warning("%s: X holds integer counts; normalized to %g and log1p'd here", sid, M.TARGET_SUM)
        log.info("%s: %d cells in the file, %d read, %d genes, %.2fs", sid, b.n_obs_file, b.X.shape[0], len(b.genes),
                 time.monotonic() - t)
        if b.X.shape[0]:
            blocks.append(b)
    if not blocks:
        log.error("every sample failed or none of the chosen cells was found")
        _write_provenance(output_dir / params["LABEL"], params, found, input_dir, statuses, notes, started, t0, {})
        return 1
    if any(b.normalized_here for b in blocks):
        notes.append("Some samples' X held integer counts; they were normalized (total 1e4, log1p) here.")
    common = set(blocks[0].genes)
    for b in blocks[1:]:
        common &= set(b.genes)
    genes = [g for g in blocks[0].genes if g in common]
    if len(genes) < len(set().union(*[set(b.genes) for b in blocks])):
        notes.append(f"Samples have different genes; the {len(genes)} genes in every sample were used.")
    lab_map = chosen.set_index(chosen["sample_id"] + "\t" + chosen["obs_id"])["label"]
    mats, obs_rows = [], []
    for b in blocks:
        mats.append(b.X[:, pd.Index(b.genes).get_indexer(genes)])
        keys = [f"{b.sample_id}\t{o}" for o in b.obs_ids]
        obs_rows.append(pd.DataFrame({"sample_id": b.sample_id, "label": lab_map.reindex(keys).to_numpy()},
                                     index=b.obs_ids))
    obs = pd.concat(obs_rows)
    adata = ad.AnnData(X=sp.vstack(mats, format="csr"), obs=obs, var=pd.DataFrame(index=pd.Index(genes)))
    del mats, blocks
    adata.obs_names_make_unique()
    counts = adata.obs["label"].value_counts()
    small = sorted(counts[counts < M.MIN_CELLS_LABEL].index)
    if small:
        notes.append(f"Not tested (fewer than {M.MIN_CELLS_LABEL} cells): " + ", ".join(small) + ".")
        adata = adata[~adata.obs["label"].isin(small)].copy()
    adata.obs["label"] = pd.Categorical(adata.obs["label"], categories=sorted(adata.obs["label"].unique()))
    adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
    t_read = time.monotonic() - t_read
    if adata.obs["label"].cat.categories.size < 2:
        log.error("need at least two labels with >= %d cells; have %s", M.MIN_CELLS_LABEL,
                  ", ".join(f"{k} ({v})" for k, v in counts.items()))
        _write_provenance(output_dir / params["LABEL"], params, found, input_dir, statuses, notes, started, t0, {})
        return 1
    log.info("%d cells x %d genes in memory (%.0f MB), %d labels, %.1fs", adata.n_obs, adata.n_vars,
             adata.X.data.nbytes / 2 ** 20 * 2, adata.obs["label"].cat.categories.size, t_read)

    # ---- test ------------------------------------------------------------------------------------
    t_test = time.monotonic()
    markers, per_label, n_genes_tested = M.find_markers(adata, params, n_threads)
    t_test = time.monotonic() - t_test

    out = output_dir / params["LABEL"]
    out.mkdir(parents=True, exist_ok=True)
    markers.to_csv(out / "markers.csv", index=False, float_format="%.6g", lineterminator="\n")

    # ---- figures ---------------------------------------------------------------------------------
    t_fig = time.monotonic()
    figures = []
    top = M.top_genes(markers, params["N_TOP"])
    if top:
        try:
            figures.append(F.dotplot(adata, top, f"Top {params['N_TOP']} markers per {lc} · {params['LABEL']}",
                                     out / "dotplot.png").name)
            figures.append(F.heatmap(adata, top, f"Top {params['N_TOP']} markers per {lc} · {params['LABEL']}",
                                     out / "heatmap.png", params["SEED"]).name)
        except Exception as exc:   # figures are a convenience; the table stands without them
            log.exception("figures failed: %s", exc)
    t_fig = time.monotonic() - t_fig

    # ---- summary / provenance --------------------------------------------------------------------
    n_samples_present = labels.groupby("label")["sample_id"].nunique()
    used = adata.obs.groupby("label", observed=True).agg(n=("sample_id", "size"), s=("sample_id", "nunique"))
    rows = []
    for _, r in per_label.iterrows():
        lab = r["label"]
        m = markers[markers["label"] == lab]
        rows.append({"label": lab, "n_cells": int(n_per_label[lab]), "n_cells_used": int(used.at[lab, "n"]),
                     "n_samples_present": int(n_samples_present[lab]), "n_samples_checked": int(r["n_samples_checked"]),
                     "sample_filter_applied": bool(r["n_samples_checked"] > 0),
                     "n_markers": len(m), "n_pooled_markers": int(r["n_pooled"]),
                     "n_dropped_inconsistent": int(r["n_dropped_inconsistent"]),
                     "top": m.sort_values("rank")["gene"].head(params["N_TOP"]).tolist()})
        log.info("label %s: %d cells (%d used) in %d samples; %d markers pass pooled, %d dropped as inconsistent "
                 "across samples, %d written", lab, n_per_label[lab], used.at[lab, "n"], n_samples_present[lab],
                 r["n_pooled"], r["n_dropped_inconsistent"], len(m))
    unfiltered = [r["label"] for r in rows if r["n_samples_checked"] == 0]
    if unfiltered:
        notes.append("Sample-aware filter not applied (no sample has >= 10 cells of the label and of the rest): "
                     + ", ".join(unfiltered) + ".")
    one_sample = [r["label"] for r in rows if r["n_samples_checked"] == 1]
    if one_sample:
        notes.append("Checked in one sample only (markers cannot be told apart from a sample effect): "
                     + ", ".join(one_sample) + ".")
    no_markers = [r["label"] for r in rows if r["n_markers"] == 0]
    if no_markers:
        notes.append("No marker passed the filters: " + ", ".join(no_markers) + ".")
    summary = {
        "label": params["LABEL"], "label_column": lc, "method": "wilcoxon",
        "n_samples": int(adata.obs["sample_id"].nunique()), "n_cells_labelled": int(len(labels)),
        "n_cells_used": int(adata.n_obs), "n_genes_common": len(genes), "n_genes_tested": n_genes_tested, "n_labels": len(rows), "n_markers": int(len(markers)),
        "n_dropped_inconsistent": int(per_label["n_dropped_inconsistent"].sum()),
        "thresholds": {"pval_adj": M.P_ADJ, "min_logfc": params["MIN_LOGFC"], "min_pct": params["MIN_PCT"],
                       "min_sample_fraction": params["MIN_SAMPLE_FRACTION"], "min_cells_sample": M.MIN_CELLS_SAMPLE},
        "labels": rows, "skipped_labels": small, "notes": notes, "figures": figures,
        "caveat": "Cells are treated as independent observations, so p-values are far too small (pseudoreplication; "
                  "Squair et al. 2021, Zimmerman et al. 2021) and, when the labels came from clustering these same "
                  "cells, circular (Zhang et al. 2019). Use them to rank genes, not as evidence; compare conditions "
                  "with st-compare's pseudobulk test.",
    }
    write_json(out / "summary.json", summary)
    timings = {"read": round(t_read, 2), "test": round(t_test, 2), "figures": round(t_fig, 2)}
    _write_provenance(out, params, found, input_dir, statuses, notes, started, t0, timings)
    log.info("done: %d labels, %d markers (%d dropped as inconsistent across samples), %.1fs, peak RSS %s MB -> %s",
             len(rows), len(markers), summary["n_dropped_inconsistent"], time.monotonic() - t0, peak_rss_mb(), out)
    for n in notes:
        log.info("note: %s", n)
    return 0


def _write_provenance(out: Path, params, found, input_dir, statuses, notes, started, t0, timings, threads=None):
    """provenance.json in the conventions' shape (plus method, skipped inputs, timings, notes, outputs)."""
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "provenance.json", {
        "app": APP_NAME, "version": app_version(), "tools": tool_versions(), "parameters": params,
        "inputs": [found.dir.rel(p) for p in found.h5ad.values()] + ([found.dir.rel(found.labels_csv)]
                                                                     if found.labels_csv else []),
        "started": started.isoformat(timespec="seconds"),
        "finished": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "samples": list(statuses.values()), "peak_rss_mb": peak_rss_mb(), "cgroup_peak_mb": cgroup_peak_mb(),
        "threads": int(os.environ.get("SC_MARKERS_THREADS") or params["THREADS"] or 1),
        "citations": CITATIONS,
        "method": {
            "cells": "per label at most MAX_CELLS_PER_LABEL cells, uniformly at random; genes present in every sample",
            "test": "scanpy.tl.rank_genes_groups(method='wilcoxon', tie_correct=True, pts=True) on log-normalized X, "
                    "each label vs all other cells, Benjamini-Hochberg per label",
            "filter": "pval_adj < 0.05, logfoldchange >= MIN_LOGFC, pct_in >= MIN_PCT, and logfoldchange >= MIN_LOGFC "
                      "within >= MIN_SAMPLE_FRACTION of the samples with >= 10 cells of the label and of the rest "
                      "(not applied to a label no sample has 10 cells of)",
            "rank": "Wilcoxon z-score among the markers that pass",
        },
        "skipped_inputs": found.skipped,
        "timings_seconds": dict(timings, total=round(time.monotonic() - t0, 2)), "notes": notes,
        "outputs": sorted(p.name for p in out.iterdir() if p.name != "provenance.json") + ["provenance.json"],
    })


if __name__ == "__main__":
    sys.exit(main())
