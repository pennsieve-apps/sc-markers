# sc-markers

Which genes mark each cell type, or each cluster? This step of the **Single-nucleus Atlas
Mapping & Condition Comparison** workflow ranks, for every label in `labels.csv`, the genes that
are higher in that label's cells than in all other cells, and keeps the ones that hold up across
samples.

`workflow A: sc-prep -> sc-qc -> sc-integrate -> sc-annotate -> {sc-markers, st-compare} -> sc-report`

**For a biologist.** For each cell type you get a ranked gene list ("NEFL, NEFM, NEFH: neurons";
"MPZ, S100B: glia"), a dot plot and a heatmap. Use it to check that a label is what its name says,
to name a cluster, or to pick genes for staining. Two cautions are built into the outputs. First,
a gene that is high in one donor's cells only is not kept as a marker: with nuclei pooled from
many donors, one donor's batch effect or contamination can otherwise top the list. Second, the
p-values are not evidence. Every nucleus counts as an independent observation although nuclei of
one donor are not independent (pseudoreplication), so p-values are far too small: Squair et al.
(2021) show the false discoveries this causes and favour pseudobulk tests; Zimmerman et al. (2021)
show the same inflation, find pseudobulk aggregation conservative and propose mixed models with a
donor random effect. When the labels are clusters of these same cells, the test also finds the
differences the clustering was built on (Zhang et al. 2019).
Read the rank, the fold change and the detection rates; for differences between conditions, use
st-compare's pseudobulk test.

**For a computational analyst.** scanpy `rank_genes_groups(method="wilcoxon", tie_correct=True,
pts=True)` is run on log-normalized `X`, each label against all other cells, pooled over samples,
with Benjamini-Hochberg per label. The Wilcoxon rank-sum test is among the most effective marker
selection methods in Pullin & McCarthy's 2024 benchmark. A **sample-aware filter** follows: in
each sample with at least 10 cells of the label and 10 of the rest, the label's log2 fold change
over the rest of *that sample* is computed (as scanpy does, from means of `expm1(X)`), and a
marker must reach `MIN_LOGFC` in at least `MIN_SAMPLE_FRACTION` of those samples. Memory is
bounded by subsampling: each label keeps at most `MAX_CELLS_PER_LABEL` cells, chosen uniformly at
random, and only those rows of each `.h5ad` are read, a block of rows at a time. That keeps the
step Lambda-sized.

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `cells` | input (`placement: subfolder`, required) | `<sample_id>.h5ad` | From sc-annotate (or sc-integrate). Only `obs_names`, `var_names` and `X` are read. `X` is log-normalized from sc-qc on; integer counts (sc-prep's output) are normalized here (total 1e4, log1p), with a note. |
| `labels` | input (`placement: subfolder`, required) | `labels.csv` (`sample_id`, `obs_id`, label columns) | `cell_type` from sc-annotate, `cluster` from sc-integrate; `obs_id` matches the h5ad's `obs_names`. |
| `markers` | output | `<LABEL>/...` | Everything under a folder named after `LABEL`, so several runs (cell types, clusters) can be linked into one port of sc-report. |

Linked, the inputs arrive in `INPUT_DIR/cells/` and `INPUT_DIR/labels/`. Without links, files
are found by name anywhere under `INPUT_DIR`. Every path comes from the shared
`processor/safe_inputs.py` (vendored unchanged from `workflow-single-cell/shared/`): symlinks are
followed only to where the platform's own links point (the entries of `INPUT_DIR` and of a real
port folder, one hop), deeper links that leave are skipped and logged; `labels.csv` is read
capped at 256 MB; each `.h5ad` is checked with `check_h5` (no external links, external storage
or virtual datasets) before it is opened. Only `.h5ad` files whose stem is a `sample_id` in
`labels.csv` are read; a `labels.csv` listing a cell twice stops the run (exit 1).

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `LABEL` | string | *(empty)* | Output folder name (`[a-z0-9-]+`). Empty: `LABEL_COLUMN` with `_` as `-` (`cell-type`, `cluster`). |
| `LABEL_COLUMN` | string | `cell_type` | Column of `labels.csv` to find markers for (`cluster` for sc-integrate's clusters). |
| `N_MARKERS` | integer | `50` | Markers written per label, best first. |
| `N_TOP` | integer | `5` | Markers per label in the dot plot and heatmap (a gene is shown once, under the first label that ranks it). |
| `MIN_PCT` | number | `0.1` | A marker is detected in at least this fraction of the label's cells. Genes below it in every label are not tested. |
| `MIN_LOGFC` | number | `0.25` | Minimum log2 fold change, pooled and within each sample that counts as agreeing (Seurat's classic `logfc.threshold`). |
| `MIN_SAMPLE_FRACTION` | number | `0.5` | Fraction of the checked samples in which the marker must reach `MIN_LOGFC`. `0` turns the filter off. |
| `MAX_CELLS_PER_LABEL` | integer | `2000` | Larger labels are subsampled to this many cells. This bounds memory and time; detection rates and fold changes of large labels are then estimated from the subsample. |
| `SEED` | integer | `42` | Seed for the subsampling and the heatmap's cell choice. |
| `THREADS` | integer | `0` | Threads; `0` = from the cgroup CPU quota, 2 without one (Lambda). |

A label needs at least 10 cells to be tested (fewer: listed in `summary.json` `skipped_labels`),
and at least two labels must remain. `p_adj < 0.05` is fixed: it is a floor on a pooled test,
not a decision.

Exit codes: `0` success; `1` no usable input (no `labels.csv` or `.h5ad`, no labelled cell in a
sample with an `.h5ad`, every sample failed, fewer than two testable labels); `2` bad parameters
(`LABEL_COLUMN` not a column, a value out of range). A sample whose file fails to read is skipped
and logged.

## Outputs (`OUTPUT_DIR/<LABEL>/`)

`markers.csv`: at most `N_MARKERS` rows per label, sorted by label and rank.

| Column | |
|---|---|
| `label`, `gene` | |
| `logfoldchange` | log2 fold change of the label over all other cells (scanpy, from means of `expm1(X)`) |
| `pval_adj` | Wilcoxon rank-sum p-value, Benjamini-Hochberg within the label. Inflated: see above. |
| `pct_in`, `pct_out` | fraction of the label's cells / of the other cells with the gene detected |
| `rank` | 1 = best, by the Wilcoxon z-score among the markers that pass every filter |
| `score` | the Wilcoxon z-score |
| `n_samples_up`, `n_samples` | samples where the label's own fold change reaches `MIN_LOGFC`; samples checked (>= 10 cells of the label and of the rest) |

The first seven columns are the contract the brief set; the last three show the ranking and the
sample-aware filter.

`dotplot.png`: the top `N_TOP` markers per label; dot size = fraction of cells, colour = mean
expression scaled per gene. `heatmap.png`: the same genes in single cells (at most 100 per label,
at random), scaled per gene.

`summary.json`: label, label column, samples, cells labelled and used, genes in common and
tested, per label (cells, cells used, samples present, samples checked, `sample_filter_applied`,
markers written, pooled hits, hits dropped as inconsistent, top `N_TOP`), skipped labels,
thresholds, notes (including "Sample-aware filter not applied" for labels no sample has 10 cells
of), figures and the p-value caveat. `provenance.json`, in the conventions' shape: `app`,
`version`, `tools`, `parameters`, `inputs` (relative paths), `started`, `finished`, `samples`
(per-sample status), `peak_rss_mb`, `cgroup_peak_mb`, `threads`, `citations` (`{doi, pmid,
text}` short form), plus `method`, `skipped_inputs`, `timings_seconds`, `notes`, `outputs`.

## Runtime

`computeTypes: [standard, lambda]`, 2 vCPU / 3008 MB, 15 min timeout. One image serves both:
`processor/entry.py` fixes the thread count (`cpu_threads`) and runs `main()` as a container, or,
when `AWS_LAMBDA_RUNTIME_API` is set, starts the Lambda runtime interface client with
`processor/handler.py`. Measured with `docker run --platform linux/amd64 --memory 3g --cpus 2
--network none --read-only --tmpfs /tmp` (**amd64 emulation on an M-series Mac**), inputs laid out
as the platform does (port folders of symlinks into the upstream output):

| Data | labels | cells used / genes tested | time (process) | peak RSS | cgroup peak |
|---|---|---|---|---|---|
| drg-c2, sc-annotate output (12 samples, 18,438 nuclei), `LABEL_COLUMN=cell_type` | 10 | 11,744 / 8,774 | 26 s | 1,209 MB | 1,139 MB |
| drg-c2, sc-annotate output, `LABEL_COLUMN=cluster` | 15 | 17,279 / 9,161 | 31 s | 1,314 MB | 1,249 MB |

Memory grows with labels x `MAX_CELLS_PER_LABEL` x genes detected, not with the dataset. An
estimate, extrapolated from the cluster run (15 labels, 17,279 cells, 1.31 GB): 30 labels at the
default cap (60,000 cells) would need about 2.5–3 GB; for more labels lower `MAX_CELLS_PER_LABEL`
or run on `standard` with more memory. Image: 929 MB unpacked (222 MB compressed; scanpy brings
numba, llvmlite, scikit-learn, statsmodels).

## Results on the demo data

drg-c2 (Discover 475; 12 samples, 18,438 nuclei after sc-qc, sc-integrate's 15 clusters,
sc-annotate's 10 cell types; acute and chronic neck pain), the upstream outputs as they stood on
2026-10-31 18:48 UTC. Outputs in the workflow repository under
`data/work/drg-c2/markers/markers/{cell-type,cluster}/`. Top five markers per cell type:

- Neuron: MAP1B, STMN2, TUBB2A, TUBB2B, SNAP25; Nociceptor (NP MRGPRD-type): PDE4C, PTGDR, KCNT1,
  ASIC3, PDIA2; Pruritogen receptor (SST+): RGS4, IL31RA, HTR3A, RET, PCSK1;
- Satellite glia: TTYH1, PRIMA1, LPL, GRIK3, LGI4; Myelinating Schwann cell: MLIP, SLC25A48,
  DRP2, KIF19, PRX;
- Fibroblast: COL1A1, COL6A2, PKDCC, DCN, COL3A1; Endothelial: VWF, ADCY4, EMCN, ECSCR, ARHGEF15;
  Macrophage: CD163, MS4A6A, FGD2, STAB1, MS4A7; Pericyte / smooth muscle: MYH11, ACTA2, CASQ2,
  CACNA1H, NOTCH3;
- Unassigned (6,463 nuclei, 35%): FOXD4L5, PPIAL4G, USP17L1, USP17L15, USP17L3, i.e. repetitive
  gene families rather than a cell type. Unassigned spans several of sc-integrate's clusters
  (C03, C04, C05 and others, each with a different best guess in sc-annotate's `clusters.csv`);
  its top markers are those of C03, the repetitive-gene cluster: likely low-quality or
  ambient-dominated nuclei that sc-qc kept (ambient removal was skipped: no GPU).

Every label got 50 markers. The sample-aware filter dropped 6,777 pooled hits for cell types (39
for Unassigned to 1,363 for Fibroblast) and 12,871 for clusters, most where a label is
concentrated in a few samples (Fibroblast in 10 samples, checked in 5; Macrophage in 8, checked
in 5). With clusters, C05 (NEFM, NEFH, NEFL with MT-CO3 and APOE) looks like a neuron-plus-ambient
mixture.

## Limitations

- **p-values are not evidence** (pseudoreplication; circularity when labels are clusters of the
  same cells). They are a floor on the list, not a test of the biology.
- **One versus rest.** A marker separates a label from everything else, not from its nearest
  neighbour. Two related subtypes share most markers, and the top genes of a small label can be
  ones it shares with a related larger type.
- **Pooled across conditions.** Markers are found on all samples together; a gene that marks a
  type only in one condition ranks lower. Use st-compare for condition effects.
- **Nuclei.** snRNA-seq measures nuclear transcripts. Some classic protein or cytoplasmic markers
  are weak, and ambient RNA from abundant types (neuronal transcripts in DRG) can appear as weak
  "markers" of other types when ambient removal did not run.
- Subsampling (`MAX_CELLS_PER_LABEL`) makes the "rest" of each label weighted towards smaller
  labels compared with their natural abundance, which changes `pct_out` somewhat.

**Best practices.** The online book that accompanies Heumos et al. (2023), *Single-cell best
practices* (sc-best-practices.org), says in its annotation chapter that "for marker gene
calculation of clusters, simple methods such as the Wilcoxon rank-sum test are thought to
perform best" (citing Pullin & McCarthy), and its differential gene expression chapter
recommends pseudobulk methods (or mixed models) for differences between conditions (checked on
2026-10-31). This app follows both. It adds the sample-aware filter, which the book does not
prescribe, because a pooled test cannot see donor effects.

## Running it locally

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python   # pytest on the synthetic fixture, three input layouts
make run PYTHON=.venv/bin/python    # builds the image (amd64), runs it over the fixture into data/output/
make validate PYTHON=.venv/bin/python   # app.yml against the manifest schema (+ CITATION.cff with cffconvert)
```

The tests check the planted markers and the dropped sample artefact in the linked, flat and
platform layouts (real port folders of symlinks into a sibling upstream output); clusters with
subsampling; counts in `X`; a missing sample and an unreadable file; a deeper symlink leaving
`INPUT_DIR`; an `.h5ad` with an HDF5 external link (refused); duplicated cells in `labels.csv`;
the "filter not applied" note; the provenance shape; bad parameters; missing inputs; a single
label; and the Lambda handler. `scripts/check-shared.sh` in the workflow repository checks that
`processor/safe_inputs.py` is the shared copy.

## Licences

Our code: Apache-2.0 (`LICENSE`). Everything in the image is under permissive licences (BSD, MIT,
Apache-2.0, PSF, MPL-2.0 unmodified); no GPL Python packages (the Debian base layer contains GPL
system tools). See `THIRD_PARTY.md`.

## Credits

- **scanpy** (Wolf, Angerer, Theis; scverse) runs the Wilcoxon test and draws the dot plot and
  heatmap; **anndata** (Virshup, Rybakov, Theis, Angerer, Wolf) holds the cells.
- The **Wilcoxon rank-sum test** (Wilcoxon 1945; Mann & Whitney 1947) and **Benjamini-Hochberg**
  (1995) are the test and the multiple-testing correction.
- **Pullin & McCarthy** (2024) benchmarked marker selection methods; their results support the
  Wilcoxon default.
- **Squair et al.** (2021), **Zimmerman et al.** (2021) and **Zhang, Kamath & Tse** (2019) showed
  why per-cell p-values are inflated or circular; their work is the caveat on every output.
- **Heumos et al.** (2023) and the Single-cell Best Practices Consortium wrote the guide this app
  follows; **Matplotlib** (Hunter 2007) draws the figures.
- Demo data: **Arendt-Tranholm A, Sankaranarayanan I and Price TJ** shared the drg-c2 single-nucleus data on Pennsieve Discover (dataset 475, CC BY 4.0);
  it is described in Arendt-Tranholm et al., *Brain* 2026.

## Citations

Each was fetched by `scripts/verify_citations.py` on 2026-10-31. Papers were checked against
Crossref (title, full author list in order, journal, year, volume, issue, pages) and PubMed
(esearch by DOI, then esummary: title, first author and year must match; PMID taken from there).
The dataset was checked against DataCite (title, creators, year, licence) and the Discover API
(DOI, contributors with ORCID, version, licence). Full author lists are in `CITATION.cff`.

- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis.
  *Genome Biol* 19(1):15 (2018). doi:10.1186/s13059-017-1382-0. PMID 29409532.
- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data
  matrices. *J Open Source Softw* 9(101):4371 (2024). doi:10.21105/joss.04371. (Crossref; not in PubMed)
- Wilcoxon F. Individual Comparisons by Ranking Methods. *Biometrics Bulletin* 1(6):80 (1945).
  doi:10.2307/3001968. (Crossref; not in PubMed)
- Mann HB, Whitney DR. On a Test of Whether one of Two Random Variables is Stochastically Larger
  than the Other. *Ann Math Stat* 18(1):50-60 (1947). doi:10.1214/aoms/1177730491. (Crossref)
- Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful
  Approach to Multiple Testing. *J R Stat Soc Series B* 57(1):289-300 (1995).
  doi:10.1111/j.2517-6161.1995.tb02031.x. (Crossref)
- Pullin JM, McCarthy DJ. A comparison of marker gene selection methods for single-cell RNA
  sequencing data. *Genome Biol* 25(1):56 (2024). doi:10.1186/s13059-024-03183-0. PMID 38409056.
- Squair JW, Gautier M, Kathe C, et al. (14 authors). Confronting false discoveries in
  single-cell differential expression. *Nat Commun* 12(1):5692 (2021).
  doi:10.1038/s41467-021-25960-2. PMID 34584091.
- Zimmerman KD, Espeland MA, Langefeld CD. A practical solution to pseudoreplication bias in
  single-cell studies. *Nat Commun* 12(1):738 (2021). doi:10.1038/s41467-021-21038-1. PMID 33531494.
- Zhang JM, Kamath GM, Tse DN. Valid Post-clustering Differential Analysis for Single-Cell
  RNA-Seq. *Cell Syst* 9(4):383-392.e6 (2019). doi:10.1016/j.cels.2019.07.012. PMID 31521605.
- Heumos L, Schaar AC, Lance C, et al. (43 authors with the consortium). Best practices for
  single-cell analysis across modalities. *Nat Rev Genet* 24(8):550-572 (2023).
  doi:10.1038/s41576-023-00586-w. PMID 37002403.
- Hunter JD. Matplotlib: A 2D Graphics Environment. *Comput Sci Eng* 9(3):90-95 (2007).
  doi:10.1109/MCSE.2007.55. (Crossref)

Demo data: Arendt-Tranholm A, Sankaranarayanan I, Price TJ. Single nucleus sequencing of human
C2 dorsal root ganglia tissue from fusion surgeries. SPARC Portal / Pennsieve Discover dataset
475, version 1 (2025). doi:10.26275/cxuf-ftav. CC BY 4.0. Described in Arendt-Tranholm A, et al.
(18 authors). Transcriptome of the human C2 dorsal root ganglia in C1–2 arthrodesis surgery:
insight for neck pain. *Brain* 149(4):1349-1364 (2026). doi:10.1093/brain/awaf370. PMID 41032656.
