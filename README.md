# FCGR-SD

**Alignment-free genome comparison on Grassmann manifolds** — a Python implementation of

> Li X., Zhou T., Feng X., Yau S.-T., Yau S. S.-T. (2024).
> *Exploring geometry of genome space via Grassmann manifolds.*
> **The Innovation** 5(5): 100677. <https://doi.org/10.1016/j.xinn.2024.100677>

English | [简体中文](README.zh-CN.md)

---

## Method

```
DNA sequence ──CGR──▶ points in the unit square ──grid counts──▶ FCGR matrix M ∈ R^{N×N}
                                                                      │
                                                                      │ column space
                                                                      ▼
                                                     a point [A] on Gr(k, N)
                                                                      │
                                                                      │ generalized Grassmann distance
                                                                      ▼
                                                      d([A₁],[A₂]) ──NJ──▶ phylogenetic tree
```

1. **CGR** maps each base to a corner of the unit square and moves the current
   point halfway toward it.
2. **FCGR** counts how many CGR points fall into each cell of a `2^n × 2^n`
   grid.
3. **Truncated SVD** keeps the leading left singular vectors, i.e. an
   orthonormal basis `A` of the column space.
4. **Principal angles** between two column spaces are computed from the SVD of
   `A₁ᵀA₂`.
5. **Generalized Grassmann distance** is the ℓ₂ norm of the principal-angle
   vector (valid across different subspace dimensions).
6. **Neighbor-joining** turns the distance matrix into a tree.

No alignment, no equal-length requirement, no k-mer frequency vector
comparison: only SVD and a few matrix products.

## Highlights

- **Exact integer grid indexing** — the CGR cell of every k-mer is computed in
  closed form, so there are no floating-point ties at cell boundaries.
- **Vectorized FCGR and truncated SVD** — one SVD per sequence; only the
  truncated basis is kept.
- **Hybrid `arccos`/`arcsin` principal angles** — accurate for nearly identical
  subspaces.
- **Complete toolkit** — NJ, Newick I/O, Robinson-Foulds distance, classical
  MDS, p-scan, and a command-line interface.
- **113 unit tests** with independent oracles (`scipy`, Biopython, brute-force
  enumeration, exact rational arithmetic).

## Installation

```bash
git clone https://github.com/ZengZichao/FCGR-SD-Python.git
cd FCGR-SD-Python
python -m venv .venv && source .venv/bin/activate
pip install -e .            # core: numpy only
# or: pip install -r requirements.txt   # includes optional test/plot/bench deps
```

Python ≥ 3.9. Core dependency: `numpy`. Optional: `scipy` (tests), `matplotlib`
(figures), `biopython` (data download and alignment baseline), `pytest`.

## Quick start

```bash
python -m fcgrsd run \
    --fasta examples/data/tableS1_ha30.fasta \
    --k 6 --p 0.71 \
    --out-prefix out/ha30
```

Outputs:

- `out/ha30.dist.csv` — distance matrix with a header and row labels;
- `out/ha30.nwk` — neighbor-joining tree in Newick format;
- `out/ha30.json` — run summary (parameters, ranks, distance statistics, timings).

## Command line

| Command | Purpose |
|---|---|
| `fcgrsd run` | FASTA → distance matrix → NJ tree |
| `fcgrsd scan-p` | Scan the truncation threshold `p` (optionally against a reference tree) |
| `fcgrsd tree` | Build an NJ tree from an existing distance matrix |
| `fcgrsd mds` | Classical MDS embedding of a distance matrix |

Common options: `--k/--n`, `--p`, `--alphabet 1..6`, `--rule ge|le`,
`--mode exact|float`, `--n-jobs`, `--clamp-negative`, `--quiet`.
Run `python -m fcgrsd <command> --help` for details.

## Python API

```python
from fcgrsd.io import read_fasta
from fcgrsd.pipeline import run_fcgr_sd

names, seqs = read_fasta("examples/data/tableS1_ha30.fasta")
res = run_fcgr_sd(names, seqs, n=6, p=0.71)

print(res.distances.shape)      # (30, 30)
print(res.ranks.min(), res.ranks.max())
print(res.tree.to_newick())
```

Lower-level building blocks:

```python
from fcgrsd.fcgr import fcgr_matrix
from fcgrsd.grassmann import subspace_basis, grassmann_distance

M1 = fcgr_matrix(seqs[0], n=6)
M2 = fcgr_matrix(seqs[1], n=6)
A1 = subspace_basis(M1, p=0.71)
A2 = subspace_basis(M2, p=0.71)
print(grassmann_distance(A1, A2))
```

## Reproducing the 30 HA example

The repository ships the 30 influenza A HA genes used in the paper's Figure 1B
(NCBI accessions from Table S1).

```bash
python -m fcgrsd run --fasta examples/data/tableS1_ha30.fasta \
    --k 6 --p 0.71 --out-prefix out/ha30
```

Expected results on this dataset:

| Quantity | Value |
|---|---|
| Optimal truncation threshold | `p = 0.71` (grid scan, RF criterion) |
| Mean truncation rank | ≈ 25 (range 24–25) |
| Monophyletic subtypes | 6 / 6 |
| Runtime (30 sequences) | ≈ 0.04 s |
| RF vs. NW+NJ reference | 4 |

The paper reports RF = 8 against a MAFFT + maximum-likelihood reference; the
reference tree here is built with Needleman–Wunsch + NJ, so the absolute RF is
not directly comparable. The qualitative result — six monophyletic subtypes —
is reproduced.

## Repository layout

```
FCGR-SD-Python/
├── fcgrsd/                 Python package
│   ├── cgr.py              CGR coordinates and exact grid indices
│   ├── fcgr.py             FCGR construction
│   ├── grassmann.py        Truncated SVD, principal angles, distance
│   ├── distance.py         Pairwise distance matrices
│   ├── tree.py             Neighbor-joining and Newick I/O
│   ├── metrics.py          Robinson-Foulds distance and classical MDS
│   ├── io.py               FASTA and distance-matrix I/O
│   ├── pipeline.py         End-to-end workflow
│   ├── viz.py              CGR/FCGR/tree/MDS plots
│   └── cli.py              Command-line interface
├── tests/                  113 unit tests
├── scripts/                Data download, simulation, figures, benchmark
├── examples/
│   ├── data/               30 HA sequences, metadata, simulated dataset
│   ├── results/            Precomputed p-scan output
│   └── figures/            Figures used in the documentation
├── docs/                   English documentation
│   └── zh-CN/              Chinese documentation
├── pyproject.toml
├── requirements.txt
├── CITATION.cff
└── LICENSE
```

## Documentation

| Document | Contents |
|---|---|
| [01 — Algorithm](docs/01-algorithm.md) | CGR/FCGR definitions, Grassmann geometry, principal angles, distance, complexity |
| [02 — Implementation](docs/02-implementation.md) | Package design, exact grid derivation, numerical details, error handling |
| [03 — Usage](docs/03-usage.md) | Installation, data preparation, CLI and API reference, parameter tuning, FAQ |
| [04 — Validation](docs/04-validation.md) | Independent checks, unit tests, reproduction results, performance, limitations |

Chinese versions: [`docs/zh-CN/`](docs/zh-CN/).

## Validation

The core is checked against independent references rather than against any
single implementation:

| Component | Reference | Result |
|---|---|---|
| Grid indexing | Exact rational CGR (`Fraction`) | identical |
| Principal angles | `scipy.linalg.subspace_angles` | max diff 9.7e-17 |
| Distance | ℓ₂ norm of SciPy angles | diff < 1e-12 |
| Truncation rank | brute-force enumeration | 0 mismatches in 1000 spectra |
| NJ | Biopython `DistanceTreeConstructor.nj` | 0 topology mismatches (100 random matrices) |
| RF | brute-force split sets | 0 mismatches (100 trees) |
| MDS | Euclidean reconstruction | max error 6.2e-15 |

Run the test suite with `python -m pytest tests -q`.

## Limitations

- **The cross-dimensional distance is not a metric.** Taking only
  `min(k₁,k₂)` principal angles can violate the triangle inequality
  (`d(A,C) > d(A,B) + d(B,C)` is possible). NJ only needs a dissimilarity, so
  tree building is fine, but branch lengths should not be interpreted as
  metric distances.
- **`p` is data-dependent.** The paper selects `p` by minimizing RF against a
  reference tree; doing this on the same data used for evaluation gives an
  optimistic estimate. Prefer a held-out subset or a fixed `p` for
  confirmatory analysis.
- **Ambiguous bases are removed.** This is the right default for draft genomes,
  but removing long `N` runs concatenates the flanking sequence and can create
  k-mers that do not exist in the genome. Split at `N` runs if this matters.
- **Negative branch lengths** can appear when the distance matrix is not
  additive. They are kept by default; use `--clamp-negative` if downstream
  tools require non-negative lengths.
- **`scan-p` caches the full left singular matrix** to avoid 99 repeated SVDs.
  That is `O(m·4^n)` memory (≈5.4 GB for 10,348 sequences at `n=8`). It is
  intended for the paper-sized datasets (30–44 sequences). The main `run`
  pipeline is streaming and keeps only the truncated bases.

## Data

- `examples/data/tableS1_ha30.fasta` — 30 influenza A HA genes downloaded from
  NCBI (public domain; accessions listed in the paper's Table S1).
- `examples/data/demo_ha30.*` — simulated sequences with a known tree,
  generated by `scripts/make_demo_data.py`.
- GISAID data are **not** redistributed. `scripts/fetch_datasets.py` can
  download sequences from NCBI given your own accession table; GISAID
  accessions are reported as missing.

## Citation

If you use this software, please cite the paper:

```bibtex
@article{Li2024FCGRSD,
  title   = {Exploring geometry of genome space via Grassmann manifolds},
  author  = {Li, X. and Zhou, T. and Feng, X. and Yau, S.-T. and Yau, S. S.-T.},
  journal = {The Innovation},
  volume  = {5},
  number  = {5},
  pages   = {100677},
  year    = {2024},
  doi     = {10.1016/j.xinn.2024.100677}
}
```

## License

MIT — see [LICENSE](LICENSE).
