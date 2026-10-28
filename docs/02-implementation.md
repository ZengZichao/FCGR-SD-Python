# Implementation

Notes on how the package is organized and why the numerical details are the way
they are. For the mathematics, see [01 — Algorithm](01-algorithm.md).

[简体中文](zh-CN/02-implementation.md) | English

---

## 1. Package layout

| Module | Responsibility |
|---|---|
| `fcgrsd.cgr` | CGR coordinates, sequence cleaning, exact/float grid indices |
| `fcgrsd.fcgr` | FCGR construction, normalization, ambiguous-base warning |
| `fcgrsd.grassmann` | Truncated SVD, principal angles, generalized distance |
| `fcgrsd.distance` | Pairwise distance matrices, optional threading, SVD cache |
| `fcgrsd.tree` | Neighbor-joining, Newick read/write, tree queries |
| `fcgrsd.metrics` | Robinson-Foulds distance, classical MDS |
| `fcgrsd.io` | FASTA and distance-matrix I/O |
| `fcgrsd.pipeline` | End-to-end workflow and p-scan |
| `fcgrsd.viz` | CGR, FCGR, tree and MDS plots |
| `fcgrsd.cli` | `python -m fcgrsd` entry point |

The dependency graph is acyclic: `cgr → fcgr → grassmann → distance → pipeline`,
with `tree`/`metrics`/`io`/`viz` used by the pipeline and CLI.

## 2. Data flow

```
FASTA ──read_fasta──▶ names, sequences
                         │
                         ▼
              fcgr_matrix(seq, n, alphabet)
                         │  (int64 count matrix, 2^n × 2^n)
                         ▼
              subspace_basis(M, p, rule)      # one SVD per sequence
                         │  (N × k orthonormal basis)
                         ▼
        distance_matrix_from_bases(bases)     # m(m-1)/2 small SVDs
                         │
                         ▼
              neighbor_joining(names, D) ──▶ Tree ──▶ Newick
```

`run_fcgr_sd` keeps one FCGR and one full `U` in memory at a time and stores only
the truncated bases, so peak memory is `O(4^n + m·4^n·k)` rather than
`O(m·4^n)`.

## 3. Exact grid indexing

The CGR recurrence gives a closed form for the grid cell (derivation in
[01](01-algorithm.md#23-closed-form-cell-index)):

```
cell = 1 + Σ_{k=0}^{n-1} bit_{i-n+1+k} · 2^k
```

The implementation encodes each base as two bit vectors `u`, `v` and accumulates
the window sums with `n` vectorized shifts:

```python
windows = len(seq) - n + 1
rows = np.ones(windows, dtype=np.int64)
for k in range(n):
    rows += u[k : k + windows] * (1 << k)
```

Complexity is `O(n·L)` per sequence, but the loop body is a single vector
operation, so the constant is small. The result is exact: no floating-point
comparison is involved, and a k-mer always maps to the same cell regardless of
its position in the sequence.

`cgr_grid_indices_float` implements the same mapping through the floating-point
recurrence (`ceil(2^n x)`, clipped to `[1, 2^n]`). It is used only as a numerical
cross-check; `mode="exact"` is the default everywhere.

## 4. FCGR construction

`fcgr_matrix` flattens the cell indices and counts them with a single
`np.bincount`:

```python
flat = (idx[:, 0] - 1) * size + (idx[:, 1] - 1)
counts = np.bincount(flat, minlength=size * size)
mat = counts.reshape(size, size)
```

**Counts vs. frequencies.** The paper calls the matrix a frequency matrix.
Dividing by the total is a positive scalar scaling, which leaves the column
space unchanged, so distances are identical. `normalize=True` is available for
users who want probabilities, but the pipeline uses counts.

**Ambiguous bases.** `normalize_sequence` removes everything outside `A/C/G/T`
and uppercases the rest. If more than 1% of a sequence is ambiguous, a warning
reports the exact fraction. Removing `N` is the right default for draft genomes,
but it concatenates the flanks of a long `N` run and can create k-mers that do
not exist in the genome; split at `N` runs if this matters.

**Resolution.** `suggest_n(L) = ceil(log₄ L)` keeps the number of cells
comparable to the sequence length. The function is a convenience, not a
constraint: `n` can be set explicitly.

## 5. Truncation and SVD reuse

The truncation rank is defined by the cumulative singular-value energy:

```python
cum = np.cumsum(s)
target = p * s.sum()
k = np.searchsorted(cum, target, side="left") + 1   # rule="ge" (default)
# k = np.searchsorted(cum, target, side="right")     # rule="le"
```

`ge` takes the smallest `k` that reaches the target; `le` takes the largest `k`
that does not exceed it. The two differ by at most one dimension. `ge` is the
default because it reproduces the published numbers; `le` follows the literal
wording of the paper. The returned rank is clamped to `[1, N]`.

**SVD cache.** A p-scan needs the same SVD for 99 values of `p`. `build_svd_cache`
stores `(U, s)` once per sequence and `basis_from_svd` slices it for each `p`.
This trades memory for time: the cache holds the full `U`, i.e. `O(m·4^n)`
(≈5.4 GB for 10,348 sequences at `n=8`). It is intended for the paper-sized
datasets (30–44 sequences). The main `run` path does **not** use the cache.

## 6. Principal angles

`principal_angles` computes `A₁ᵀA₂ = UΣVᵀ` and splits the angles at
`ε = cos(1e-5)`:

- angles with `σ ≥ ε` (small angles) are computed as `arcsin(μ)` from the
  projected matrix `B = (I − A₁A₁ᵀ)A₂V`;
- the remaining angles are computed as `arccos(σ)`.

Both `σ` and `μ` are clipped to `[0, 1]` before the inverse trigonometric
functions. Rounding can produce `σ = 1 + 9e-16`; without clipping, `arccos`
returns a complex number and the whole distance matrix is promoted to a complex
dtype.

The small-angle block is sorted ascending before `arcsin`, so the returned
vector satisfies `θ₁ ≤ θ₂ ≤ …`. The distance is a norm and therefore
permutation-invariant, but the ordering is part of the public contract and is
checked by a regression test.

## 7. Distance matrix

`distance_matrix_from_bases` fills the upper triangle and mirrors it:

```python
for i in range(m - 1):
    for j in range(i + 1, m):
        d = grassmann_distance(bases[i], bases[j])
        D[i, j] = D[j, i] = d
```

Each unordered pair is computed once, so the result is exactly symmetric. With
`n_jobs > 1`, rows are distributed over a `ThreadPoolExecutor`; LAPACK releases
the GIL during SVD, so this gives real parallelism without copying the bases.
Each thread writes distinct cells of `D`.

## 8. Neighbor-joining and Newick

The implementation is the classical Saitou–Nei algorithm:

```
while more than two nodes remain:
    r_i = Σ_k D_ik
    Q_ij = (n-2) D_ij - r_i - r_j
    merge the pair minimizing Q
    d(i,u) = D_ij/2 + (r_i - r_j)/(2(n-2))
    d(j,u) = D_ij - d(i,u)
    D_uk = (D_ik + D_jk - D_ij)/2
```

Distances are updated in place; dead rows/columns are zeroed and masked, which
avoids copying submatrices. The last two nodes are joined with half the remaining
distance each, giving a degree-2 root. Ties are broken by the flattened
row-major argmin, so the result is deterministic.

NJ can produce negative branch lengths when the matrix is not additive. They are
kept by default; `clamp_negative=True` (CLI `--clamp-negative`) truncates them to
zero.

Newick output quotes labels containing whitespace or `()[],:;'"` and strips
trailing zeros from branch lengths. `read_newick` handles quoted labels,
comments in `[...]`, and branch lengths.

## 9. Robinson-Foulds and MDS

`robinson_foulds` compares unrooted bipartitions by default and rooted clades
with `rooted=True`. Trivial splits (one leaf, or all but one) are ignored. The
normalized value is `|S_a Δ S_b| / (|S_a| + |S_b|)`, which lies in `[0, 1]`.
Trees with duplicate leaf names are rejected: the split sets would collapse.

`classical_mds` is Torgerson's PCoA:

```python
B = -0.5 * J @ (D**2) @ J
vals, vecs = np.linalg.eigh(B)
coords = vecs[:, :k] * np.sqrt(np.clip(vals[:k], 0, None))
```

Negative eigenvalues are clipped to zero, which is the standard PCoA treatment
for non-Euclidean dissimilarities. The number of returned components is
`min(n_components, m)`.

## 10. Error handling

| Input | Behavior |
|---|---|
| Sequence shorter than `n` | `ValueError` |
| Empty or all-ambiguous sequence | `ValueError` |
| `n` not a positive integer | `ValueError` |
| `p` outside `(0, 1)` | `ValueError` |
| `alphabet` outside `1..6` | `ValueError` |
| Duplicate sequence names | `ValueError` (FASTA reader, pipeline, NJ) |
| Distance matrix with `NaN`/`Inf` | `ValueError` |
| Distance-matrix header ≠ row labels | `ValueError` |
| Fewer than two sequences in `run` | `ValueError` |
| `mds --dims` larger than the number of points | clipped to `m` |

The goal is to fail loudly on inputs that would otherwise produce a plausible
but meaningless tree.

## 11. Memory and performance

- `run_fcgr_sd`: one FCGR + one full `U` at a time; retains only truncated bases
  (`O(m·4^n·k)`).
- `scan_truncation_p`: retains full `U` for all sequences (`O(m·4^n)`) to avoid
  99 repeated SVDs.
- Pairwise distances: `O(m²)` output; `k³` SVD per pair plus an `N·k²`
  projection.
- Threading: `--n-jobs` parallelizes the pair loop; no process pool, so bases are
  shared without serialization.

Measured on the 30 HA dataset (`n=6`, `p=0.71`): ≈0.04 s end to end, mean rank
≈25.
