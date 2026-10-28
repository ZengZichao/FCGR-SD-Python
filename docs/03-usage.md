# Usage

[简体中文](zh-CN/03-usage.md) | English

---

## 1. Installation

```bash
git clone https://github.com/ZengZichao/FCGR-SD-Python.git
cd FCGR-SD-Python
python -m venv .venv && source .venv/bin/activate
pip install -e .
```

Requirements: Python ≥ 3.9, `numpy ≥ 1.22`. Optional extras:

| Extra | Packages | Used by |
|---|---|---|
| tests | `pytest`, `scipy` | `pytest tests` |
| plots | `matplotlib` | `scripts/make_figures.py`, `fcgrsd mds --plot` |
| bench | `scipy`, `biopython` | `scripts/benchmark_alignment.py`, `scripts/fetch_datasets.py` |

Install everything with `pip install -r requirements.txt`.

Check the installation:

```bash
python -m pytest tests -q      # expected: 113 passed
python -m fcgrsd --help
```

## 2. Data preparation

### FASTA

One record per sequence, unique names, standard IUPAC characters. Ambiguous
bases (`N`, `-`, etc.) are removed before CGR; a warning is emitted when they
exceed 1% of a sequence.

```
>KM409247.1
ATGAAGGCAATACTAGTAGTTCTGCTATACACATTTGCAACCACAAATGCAGACACATTATGTATAGGTT...
```

Names are taken from the first whitespace-delimited token of the header, so
`>seq1 description` gets the name `seq1`.

### Accession tables

`scripts/fetch_datasets.py` accepts a CSV with a header row. The accession column
is detected from headers containing `accession`, `gisaid`, `header` or `id`; a
`Header` value such as `"LC710779.1 Avian coronavirus"` is reduced to
`LC710779.1`. An optional group column is detected from `gene`, `type`,
`variant`, `subgenus`, `lineage`, `pangolin` or `clade`.

```bash
python scripts/fetch_datasets.py \
    --table my_accessions.csv \
    --out data/sequences.fasta \
    --batch 50
```

GISAID accessions (`EPI_ISL_*`) are usually not available from NCBI and are
written to `*_missing.txt`.

## 3. Command line

### `run` — FASTA to tree

```bash
python -m fcgrsd run \
    --fasta examples/data/tableS1_ha30.fasta \
    --k 6 --p 0.71 \
    --out-prefix out/ha30
```

Progress output (the CLI messages are currently in Chinese; the values are the
important part):

```
[Step 1] 计算 FCGR 与截断 SVD ...
         完成，用时 0.015s；截断维数 min/median/max = 24/25/25
[Step 3] 计算成对广义 Grassmann 距离 ...
         完成 435 对，用时 0.027s
[Step 4] Neighbor-Joining 建树 ...
         完成，用时 0.001s
[完成] 总用时 0.043s
已写出：out/ha30.dist.csv / out/ha30.nwk / out/ha30.json
```

| Option | Default | Meaning |
|---|---|---|
| `--fasta` | required | input FASTA |
| `--k` / `--n` | `ceil(log₄ L)` | FCGR resolution, matrix `2^k × 2^k` |
| `--p` | 0.75 | cumulative singular-value threshold |
| `--alphabet` | 3 | base-to-corner convention, 1–6 |
| `--rule` | `ge` | truncation criterion: `ge` or `le` |
| `--mode` | `exact` | grid indexing: `exact` or `float` |
| `--n-jobs` | 1 | threads for the pairwise distance loop |
| `--clamp-negative` | off | truncate negative NJ branch lengths to 0 |
| `--out-prefix` | FASTA basename | output prefix |
| `--quiet` | off | suppress progress output |

Outputs:

- `.dist.csv` — header + row labels + 10-decimal distances;
- `.nwk` — Newick tree;
- `.json` — parameters, rank statistics, distance statistics, timings.

### `scan-p` — choose the truncation threshold

```bash
python -m fcgrsd scan-p \
    --fasta examples/data/tableS1_ha30.fasta \
    --k 6 \
    --reference examples/results/tableS1_scan_ref.nwk \
    --out out/scan.csv
```

Scans `p ∈ [0.01, 0.99]` in steps of 0.01 and writes one row per `p`:
`p, rank_mean, rank_min, rank_max, rf, rf_normalized`. With `--reference`, it
prints the `p` with the smallest RF (ties broken toward smaller `p`).

> `scan-p` caches the full `U` of every sequence to avoid 99 repeated SVDs.
> That is `O(m·4^n)` memory; it is intended for tens of sequences, not thousands.

### `tree` — distance matrix to tree

```bash
python -m fcgrsd tree --dist out/ha30.dist.csv --out out/ha30.nwk
python -m fcgrsd tree --dist out/ha30.dist.csv --out out/ha30.nwk --clamp-negative
```

Useful when the distance matrix comes from another tool.

### `mds` — classical MDS

```bash
python -m fcgrsd mds --dist out/ha30.dist.csv --dims 3 \
    --out out/mds.csv --plot out/mds.png
```

`--dims` larger than the number of points is clipped to `m`, so the CSV header
always matches the number of coordinate columns.

## 4. Python API

### End-to-end

```python
from fcgrsd.io import read_fasta
from fcgrsd.pipeline import run_fcgr_sd

names, seqs = read_fasta("examples/data/tableS1_ha30.fasta")
res = run_fcgr_sd(names, seqs, n=6, p=0.71)

res.distances          # (30, 30) symmetric matrix
res.ranks              # truncation rank per sequence
res.tree               # Tree object
res.timings            # dict of stage timings
res.summary()          # JSON-serializable summary
```

`run_fcgr_sd` accepts `rule="ge"|"le"`, `mode="exact"|"float"`,
`alphabet=1..6`, `n_jobs`, `build_tree`, `clamp_negative`.

### Building blocks

```python
from fcgrsd.fcgr import fcgr_matrix
from fcgrsd.grassmann import subspace_basis, principal_angles, grassmann_distance

M1 = fcgr_matrix(seqs[0], n=6)
M2 = fcgr_matrix(seqs[1], n=6)

A1 = subspace_basis(M1, p=0.71)          # (64, k1) orthonormal basis
A2 = subspace_basis(M2, p=0.71)

theta = principal_angles(A1, A2)         # ascending principal angles
d = grassmann_distance(A1, A2)           # ||theta||_2
```

### Trees and metrics

```python
from fcgrsd.tree import neighbor_joining, read_newick
from fcgrsd.metrics import robinson_foulds, classical_mds

tree = neighbor_joining(names, res.distances)
reference = read_newick(open("reference.nwk").read())
print(robinson_foulds(tree, reference))
print(robinson_foulds(tree, reference, normalize=True))

coords, eigvals = classical_mds(res.distances, 3)
```

### Scanning `p`

```python
from fcgrsd.pipeline import scan_truncation_p

rows = scan_truncation_p(
    names, seqs, n=6,
    ps=[round(0.01 * i, 2) for i in range(1, 100)],
    reference_tree=reference,
)
best = min(rows, key=lambda r: (r["rf"], r["p"]))
print(best)
```

## 5. Parameter tuning

### `n` — resolution

Keep the number of cells `4^n` comparable to the sequence length:

| Sequence length | Suggested `n` | Matrix |
|---|---|---|
| < 500 | 4 | 16 × 16 |
| 1–3 kb (HA, spike) | 6 | 64 × 64 |
| 10–40 kb (SARS-CoV-2) | 8 | 256 × 256 |
| > 100 kb (bacterial genomes) | 9–10 | 512²–1024² |

Too small loses k-mer detail; too large makes the matrix sparse and the SVD
slow. `n > 11` triggers a memory warning.

### `p` — truncation

`p` controls how much of the singular-value spectrum is kept. The paper's values
are 0.71 for the 30 HA genes and 0.75 for the 44 SARS-CoV-2 genomes. For a new
dataset, scan `p` against a reference tree or a held-out subset. Do not tune and
evaluate on the same data if you report the RF as a performance estimate.

### `alphabet`

`alphabet=3` (AGCT) is the paper's convention and the default. Other assignments
rotate/reflect the CGR square; they are useful as a sensitivity check but change
the FCGR row/column structure.

### `rule`

`ge` (default) takes the smallest `k` whose cumulative energy reaches `p·total`;
`le` takes the largest `k` that does not exceed it. They differ by at most one
dimension. Use `le` only if you need to follow the paper's wording literally.

### `mode`

`exact` is the default and should be used for production. `float` follows the
floating-point CGR recurrence and is useful only for numerical comparison.

## 6. FAQ

**The tree has negative branch lengths.**
NJ can produce them when the distance matrix is not additive. This is expected;
use `--clamp-negative` if a downstream tool requires non-negative lengths.

**A sequence contains many `N`s.**
They are removed, with a warning above 1%. If `N` occurs in long runs, removing
it concatenates the flanks and can create spurious k-mers; split the sequence at
`N` runs and count the pieces separately.

**Memory is exhausted.**
The main `run` path is streaming. If `scan-p` runs out of memory, use fewer
sequences or a smaller `n`; the p-scan cache is `O(m·4^n)`.

**Can I use protein sequences?**
Not directly. The method relies on the four-corner geometry of DNA; a protein
version needs a 20-vertex construction and a different grid encoding.

**The RF differs from the paper.**
The paper uses MAFFT + maximum likelihood as the reference; this package uses
Neighbor-Wunsch + NJ in the benchmark script. Absolute RF values are therefore
not comparable, but the p-selection protocol and the qualitative conclusions are.

**Why does `run` require at least two sequences?**
A single sequence has no pairwise distance and no meaningful tree. Use
`distance_matrix` directly if you only need a 1×1 matrix.

**Can I use a custom distance matrix?**
Yes: `fcgrsd tree --dist matrix.csv` and `fcgrsd mds --dist matrix.csv`. The CSV
must have a header, row labels, and a square symmetric matrix.
