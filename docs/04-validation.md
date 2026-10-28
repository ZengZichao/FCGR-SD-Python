# Validation

[简体中文](zh-CN/04-validation.md) | English

---

## 1. Strategy

Correctness is checked against **independent references**, not against a single
implementation:

- exact rational arithmetic (`fractions.Fraction`) for the CGR grid;
- `scipy.linalg.subspace_angles` for principal angles;
- brute-force enumeration for the truncation rank and Robinson-Foulds distance;
- Biopython's neighbor-joining for tree topology;
- Euclidean reconstruction for MDS.

The goal is to catch both mathematical errors and implementation slips
(off-by-one indexing, wrong ordering, silent dtype promotion).

## 2. Test suite

```bash
python -m pytest tests -q     # 113 passed
```

| File | Tests | Coverage |
|---|---|---|
| `test_cgr.py` | 42 | grid indices vs exact rational CGR (6 alphabets × 5 resolutions); float vs integer paths; degenerate homopolymers; coordinate range; invalid inputs |
| `test_fcgr.py` | 16 | count conservation; single-k-mer cell formula; normalization; ambiguous-base removal and warning; lowercase equivalence; `fcgr_total` |
| `test_grassmann.py` | 15 | principal angles vs SciPy; small-angle accuracy and ascending order; triangle inequality; rotation invariance; truncation monotonicity; parameter checks |
| `test_tree.py` | 15 | NJ recovers additive trees; Newick round-trip; known RF values; quoted names; split counts; NaN/duplicate rejection; negative-branch clamping |
| `test_pipeline.py` | 10 | end-to-end structure; identical sequences; distance monotonicity; `ge` vs `le`; MDS reconstruction; duplicate/single-sequence rejection |
| `test_io.py` | 5 | distance-matrix round-trip; header/label consistency; square check; empty and duplicate FASTA records |
| `test_scripts.py` | 6 | accession-column detection for five table layouts |
| `test_viz.py` | 2 | FCGR heatmap orientation matches CGR; zero-length tree does not produce a degenerate axis |
| `test_cli.py` | 2 | `mds --dims > m` header; `tree --clamp-negative` |

## 3. Numerical checks

### CGR grid

For six alphabets and `n ∈ {1,2,3,6,8}`, the integer closed form is compared
with a `Fraction`-based exact CGR followed by `floor(2^n x) + 1`. All cells are
identical. The floating-point path agrees with the integer path on random and
homopolymer sequences.

### Principal angles

On random subspaces (including four angles of order `1e-7`), the maximum
difference from `scipy.linalg.subspace_angles` is `9.7e-17`, and the returned
angles are ascending.

### Distance

For subspaces of dimensions 4 and 6, `grassmann_distance` equals the ℓ₂ norm of
SciPy's principal angles to within `1e-12`.

### Truncation

For 1,000 random singular-value spectra and random `p`, `truncation_rank`
matches a brute-force scan for both `ge` and `le` rules.

### Neighbor-joining

- 200 random additive trees: the cophenetic matrix is recovered to `3.6e-15`,
  and RF against the true tree is 0.
- 100 random non-additive Euclidean distance matrices: the topology matches
  Biopython's NJ in all cases.
- Serial and 4-thread distance matrices are bitwise identical.

### Robinson-Foulds

For 100 random tree pairs, the split-based RF equals a brute-force symmetric
difference of bipartitions.

### MDS

For random Euclidean configurations, classical MDS reconstructs the distances to
`6.2e-15`.

## 4. Reproduction on the 30 HA dataset

Data: `examples/data/tableS1_ha30.fasta` (30 influenza A HA genes, NCBI
accessions from the paper's Table S1).

```bash
python -m fcgrsd run --fasta examples/data/tableS1_ha30.fasta \
    --k 6 --p 0.71 --out-prefix out/ha30
```

| Quantity | Result |
|---|---|
| Optimal `p` (0.01–0.99 scan, RF criterion) | **0.71** |
| Truncation rank | 24–25 (mean ≈ 25) |
| Monophyletic subtypes | **6 / 6** |
| Runtime | **0.039 s** |
| RF vs. NW+NJ reference | **4** |

The six subtypes (H1N1, H2N2, H3N2, H5N1, H7N3, H7N9) each form a monophyletic
group, reproducing the qualitative result of the paper's Figure 1B. The paper
reports RF = 8 against a MAFFT + maximum-likelihood reference; the reference
tree here is Needleman–Wunsch + NJ, so the absolute value is not comparable.

The p-scan output is in `examples/results/tableS1_scan.csv`; the curve is shown
in `examples/figures/fig4_scan_p.png`.

## 5. Simulated data

`examples/data/demo_ha30.*` contains 30 simulated sequences (1,700 bp) generated
along a known random tree with Jukes–Cantor substitutions. Running the pipeline
with `n=6, p=0.71` gives RF = 8 against the true tree (normalized 0.148), which
is consistent with the expected loss of resolution on deeply diverged simulated
branches.

```bash
python scripts/make_demo_data.py --out-dir examples/data --seed 20240909
python -m fcgrsd run --fasta examples/data/demo_ha30.fasta \
    --k 6 --p 0.71 --out-prefix out/demo
```

## 6. Performance and memory

Measured on the 30 HA dataset (`n=6`, `p=0.71`, one core):

| Stage | Time |
|---|---|
| FCGR + truncated SVD | 0.015 s |
| 435 pairwise distances | 0.027 s |
| NJ | 0.0005 s |
| **Total** | **0.039 s** |

Memory:

- `run`: one FCGR and one full `U` at a time; retained bases are `O(m·4^n·k)`;
- `scan-p`: full `U` for every sequence, `O(m·4^n)` (≈5.4 GB for 10,348
  sequences at `n=8`), intended for tens of sequences.

## 7. Limitations

The following are properties of the method, not implementation bugs, and are
documented so that results are not over-interpreted:

1. **Cross-dimensional distance is not a metric.** Taking only `min(k₁,k₂)`
   principal angles can violate the triangle inequality. NJ only needs a
   dissimilarity; branch lengths are not metric distances.
2. **`p` selection is data-dependent.** Choosing `p` by minimizing RF on the same
   data used for evaluation gives an optimistic estimate.
3. **Distance saturation.** Highly divergent sequences have nearly orthogonal
   subspaces; distances approach `√k·π/2` and deep branches lose resolution.
4. **Ambiguous bases.** Removing `N` runs concatenates flanks and can create
   chimeric k-mers.
5. **No confidence values.** The pipeline returns a point estimate (a tree); it
   does not provide bootstrap support.

## 8. Reproducing the results

```bash
# Install
pip install -r requirements.txt

# Unit tests
python -m pytest tests -q

# 30 HA example
python -m fcgrsd run --fasta examples/data/tableS1_ha30.fasta \
    --k 6 --p 0.71 --out-prefix out/ha30

# p scan against the NW+NJ reference
python -m fcgrsd scan-p --fasta examples/data/tableS1_ha30.fasta --k 6 \
    --reference examples/results/tableS1_scan_ref.nwk --out out/scan.csv

# Figures (requires matplotlib)
python scripts/make_figures.py \
    --fasta examples/data/tableS1_ha30.fasta \
    --meta examples/data/tableS1_ha30_meta.csv \
    --k 6 --p 0.71 \
    --scan examples/results/tableS1_scan.csv \
    --out-dir examples/figures
```

> The reference tree `examples/results/tableS1_scan_ref.nwk` is the
> Needleman–Wunsch + NJ tree used for the p scan. Regenerate it with
> `scripts/benchmark_alignment.py` if you change the input data.
