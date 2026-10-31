# Third-party components

Our code (`processor/`, `tests/`, `scripts/`) is Apache-2.0 (`LICENSE`). The Docker image bundles
the components below. Licences were read from each package's installed metadata
(`License-Expression`, `License`, classifiers) in the built image, and from its licence file where the
metadata says only "BSD" (umap-learn, numba, seaborn, kiwisolver), on 2026-10-31.

| Component | Version | Licence |
|---|---|---|
| Debian (bookworm, via `python:3.12.16-slim-bookworm`) | 12 | various free licences (per package, `/usr/share/doc/*/copyright`) |
| CPython | 3.12.16 | PSF-2.0 |
| scanpy | 1.12.5 | BSD-3-Clause |
| anndata | 0.13.4 | BSD-3-Clause |
| h5py | 3.16.0 | BSD-3-Clause (bundled HDF5 under its own BSD-style licence) |
| NumPy | 2.5.4 | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0 (bundled OpenBLAS: BSD-3-Clause) |
| SciPy | 1.18.1 | BSD-3-Clause |
| pandas | 3.1.0 | BSD-3-Clause |
| Matplotlib | 3.11.2 | Matplotlib License (PSF-based); bundles DejaVu and STIX fonts under their own free licences |
| scikit-learn | 1.9.1 | BSD-3-Clause |
| statsmodels | 0.15.0 | BSD-3-Clause |
| seaborn | 0.13.2 | BSD-3-Clause |
| numba | 0.68.1 | BSD-2-Clause |
| llvmlite | 0.50.0 | BSD-2-Clause AND Apache-2.0 WITH LLVM-exception |
| umap-learn | 0.5.12 | BSD-3-Clause |
| pynndescent | 0.7.0 | BSD-2-Clause |
| networkx | 3.7 | BSD-3-Clause |
| patsy | 1.0.3 | BSD-2-Clause |
| joblib | 1.6.0 | BSD-3-Clause |
| threadpoolctl | 3.7.0 | BSD-3-Clause |
| cloudpickle | 3.1.2 | BSD-3-Clause |
| contourpy | 1.4.0 | BSD-3-Clause |
| cycler | 0.12.1 | BSD-3-Clause |
| kiwisolver | 1.5.1 | BSD-3-Clause ("Modified BSD") |
| msgspec | 0.22.0 | BSD-3-Clause |
| python-dotenv | 1.2.4 | BSD-3-Clause |
| scverse-misc | 0.1.8 | BSD-3-Clause |
| wrapt | 2.5.2 | BSD-2-Clause |
| zarr | 3.4.2 | MIT |
| numcodecs | 0.17.1 | MIT |
| fonttools | 4.68.0 | MIT |
| formulaic | 1.2.2 | MIT |
| interface_meta | 2.0.1 | MIT |
| narwhals | 2.28.0 | MIT |
| natsort | 8.4.0 | MIT |
| array-api-compat | 1.16.0 | MIT |
| donfig | 0.8.1.post1 | MIT |
| pydantic, pydantic_core, pydantic-settings | 2.14.2, 2.50.2, 2.16.2 | MIT |
| annotated-types, typing-inspection | 0.8.0, 0.4.4 | MIT |
| pyparsing | 3.3.3 | MIT |
| PyYAML | 6.0.3 | MIT |
| six | 1.17.0 | MIT |
| Pillow | 13.0.0 | MIT-CMU |
| typing_extensions | 4.16.0 | PSF-2.0 |
| fast-array-utils | 1.5.1 | MPL-2.0 |
| legacy-api-wrap | 1.5 | MPL-2.0 |
| session-info2 | 0.4.2 | MPL-2.0 |
| tqdm | 4.70.1 | MPL-2.0 AND MIT |
| certifi | 2026.10.23 | MPL-2.0 |
| packaging | 26.3 | Apache-2.0 OR BSD-2-Clause |
| python-dateutil | 3.0.0 | Apache-2.0 OR BSD-3-Clause |
| google-crc32c | 1.9.0 | Apache-2.0 |
| tzdata | 2026.6 | Apache-2.0 |
| awslambdaric | 4.2.3 | Apache-2.0 |
| snapshot-restore-py (awslambdaric dependency) | 1.0.0 | Apache-2.0 |
| simplejson (awslambdaric dependency) | 4.2.0 | MIT OR AFL-2.1 |
| pip | 25.0.1 | MIT |

No GPL Python packages are in the image (the Debian base layer contains GPL system tools, e.g.
bash and coreutils): scanpy's optional Leiden backends (igraph, leidenalg) are not
installed, as this app does no clustering. The MPL-2.0 components are used unmodified. The demo
data the README describes (Discover 475, CC BY 4.0) is not in the image.
