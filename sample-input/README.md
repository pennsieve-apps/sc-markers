# sample-input

This app reads two inputs (the contract is in `workflow-single-cell/docs/conventions.md`):

- `cells/<sample_id>.h5ad` from `sc-annotate` (or `sc-integrate`): log-normalized `X`
  (`sc-prep`'s integer counts also work: they are normalized here), `obs_names` =
  `<sample_id>:<barcode>`,
- `labels/labels.csv` from `sc-annotate` (`cell_type`) or `sc-integrate` (`cluster`).

## Synthetic fixture (no download)

```sh
make fixture PYTHON=.venv/bin/python   # data/fixture/{linked,flat}/
make run PYTHON=.venv/bin/python       # Docker over data/fixture/linked/ -> data/output/cell-type/
```

Four samples (`s1`, `s2` control; `s3`, `s4` pain), 400 nuclei and 300 genes each; three cell
types (Neuron 15 %, Satellite_glia 35 %, Fibroblast 50 %) with five planted markers each (4x),
clusters `C01`-`C04`. Planted artefact: `SAMPLE_GENE` is 12x higher in `s1`'s neurons only. The
app ranks every planted marker in its type's top five and drops `SAMPLE_GENE` (high in one of four
samples); with `MIN_SAMPLE_FRACTION=0` it is kept.

## Real data

The `drg-c2` demo of the workflow repository (Discover 475, CC BY 4.0): run
`sc-prep -> sc-qc -> sc-integrate -> sc-annotate`, then give this app sc-annotate's `cells` and
`labels` outputs, one folder each:

```sh
W=../workflow-single-cell/data/work/drg-c2
mkdir -p in && cp -R $W/annotate/cells in/cells && cp -R $W/annotate/labels in/labels
make run INPUT=$PWD/in
```

Mount real copies, not symlinks to host paths (they do not resolve inside the container).
