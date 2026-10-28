# Algorithm

This document gives the mathematical background of FCGR-SD, following the
paper's main text and supplementary §1.2 and §2–§4. Proofs are included where
they clarify the implementation.

[简体中文](zh-CN/01-algorithm.md) | English

---

## 1. Overview

Alignment-based comparison scales poorly with the number of sequences and their
length. FCGR-SD is **alignment-free**: each sequence is compressed into a
matrix, the matrix is mapped to a point on a manifold, and similarity is
measured by an intrinsic distance on that manifold.

```
DNA sequence ──CGR──▶ points in the unit square ──grid counts──▶ FCGR M ∈ R^{N×N}
                                                                     │
                                                                     │ column space
                                                                     ▼
                                                     a point [A] on Gr(k, N)
                                                                     │
                                                                     │ generalized Grassmann distance
                                                                     ▼
                                                      d([A₁],[A₂]) ──NJ──▶ tree
```

Here `N = 2^n`, and `n` is the FCGR resolution (called `k` in some parts of the
paper).

---

## 2. Chaos Game Representation (CGR)

### 2.1 Definition

Assign the four bases to the corners of the unit square. Starting from the
center `P₀ = (1/2, 1/2)`, move the current point halfway toward the corner of
each new base:

$$
P_0=\left(\tfrac12,\tfrac12\right),\qquad
P_i=\tfrac12\bigl(P_{i-1}+P(s_i)\bigr),\quad i=1,\dots,L .
$$

With the `AGCT` convention used by the main pipeline:

| Base | Corner `P(s)` |
|---|---|
| A | (0, 0) |
| G | (0, 1) |
| C | (1, 1) |
| T | (1, 0) |

The package also supports five other corner assignments (`alphabet = 1..6`):

| `alphabet` | Order | Corners `(x, y)` |
|---|---|---|
| 1 | ACGT | A(0,0), C(0,1), G(1,1), T(1,0) |
| 2 | ACTG | A(0,0), C(0,1), T(1,1), G(1,0) |
| 3 | AGCT | A(0,0), G(0,1), C(1,1), T(1,0) — **default; paper Algorithm 2** |
| 4 | AGTC | A(0,0), G(0,1), T(1,1), C(1,0) |
| 5 | ATCG | A(0,0), T(0,1), C(1,1), G(1,0) |
| 6 | ATGC | A(0,0), T(0,1), G(1,1), C(1,0) |

Different assignments rotate/reflect the point cloud; they do not change the
underlying k-mer statistics, but at finite resolution they change the row/column
structure of the FCGR and therefore the truncated subspace.

### 2.2 Binary structure

**Lemma 1.** Let `u_i ∈ {0,1}` be the first coordinate of the i-th base. Then

$$
x_i=\frac{X_i}{2^{\,i+1}},\qquad X_i=1+\sum_{j=1}^{i}u_j\,2^{j},
$$

and `X_i` is always **odd**.

*Proof.* By induction. `x₀ = 1/2 = X₀/2` with `X₀ = 1`. If
`x_{i-1} = X_{i-1}/2^i`, then

$$
x_i=\frac{x_{i-1}+u_i}{2}=\frac{X_{i-1}+u_i 2^{i}}{2^{\,i+1}}
\;\Longrightarrow\;
X_i=X_{i-1}+u_i2^i .
$$

`X_{i-1}` is odd and `u_i 2^i` is even for `i ≥ 1`, so `X_i` is odd. ∎

### 2.3 Closed-form cell index

Divide the unit square into `2^n × 2^n` cells numbered from 1. With
`m = i+1-n ≥ 1`:

$$
\left\lfloor 2^n x_i\right\rfloor
=\left\lfloor \frac{X_i}{2^{m}}\right\rfloor
=\left\lfloor \frac{\overbrace{1+\sum_{j<m}u_j2^j}^{<\,2^m}+\sum_{j=m}^{i}u_j2^j}{2^m}\right\rfloor
=\sum_{k=0}^{n-1}u_{\,i-n+1+k}\,2^{k}.
$$

Therefore

$$
\boxed{\;c_i^{(1)}=1+\sum_{k=0}^{n-1}u_{\,i-n+1+k}2^{k},\qquad
c_i^{(2)}=1+\sum_{k=0}^{n-1}v_{\,i-n+1+k}2^{k}\;}
$$

Three consequences:

1. The cell depends only on the **last n bases of the window** — matching the
   standard n-mer counting interpretation.
2. Since `X_i` is odd, `2^n x_i` never lies exactly on a grid line, so
   `ceil(2^n x_i) = floor(2^n x_i) + 1` in exact arithmetic.
3. The formula is **pure integer arithmetic** and avoids IEEE-754 rounding
   entirely. `mode="exact"` (default) uses it; `mode="float"` follows the
   floating-point recurrence for numerical comparison.

### 2.4 FCGR and k-mer counts

Because the four bases have distinct `(u, v)` values, the map from n-mers to
cells is a bijection `4^n → 4^n`.

> **Proposition 1.** The FCGR is a rearrangement (row/column permutation) of the
> n-mer count table. FCGR-SD is therefore still a k-mer composition method, but
> it compares the **column space** of the count matrix rather than the frequency
> vector itself.

This explains why the method is insensitive to sequence length (all sequences
share the same ambient space `R^N`) and why it is more robust than
"FCGR + Euclidean/Pearson distance": the latter is sensitive to individual cell
counts, while the former only cares about the position of a subspace.

---

## 3. FCGR

Count the points in each cell:

$$
M=(m_{ij})\in\mathbb{R}^{N\times N},\qquad
m_{ij}=\#\{i': c^{(1)}_{i'}=i,\;c^{(2)}_{i'}=j\},
$$

with `N = 2^n` and `L - n + 1` points in total. The paper calls this a
"frequency matrix"; storing raw counts is equivalent because the column space is
invariant under positive scaling.

**Resolution.** A good rule of thumb is `n ≈ log₄ L`, so that the number of
cells `4^n` is comparable to the sequence length. The paper's two examples
follow it:

- 30 HA genes, mean length 1,718.8 → `n = 6` (`4⁶ = 4096`);
- 44 SARS-CoV-2 genomes, mean length 29,759 → `n = 8` (`4⁸ = 65536`).

The package provides `fcgrsd.pipeline.suggest_n(L) = ceil(log₄ L)`.

---

## 4. Why the column space

The column space of `M ∈ R^{N×N}` is a subspace of `R^N`. Two reasons motivate
this choice:

1. **Biological**: columns correspond to the second CGR coordinate; taking the
   column space treats the distribution along a column as the feature, which is
   compatible with the purine/pyrimidine geometry of the square.
2. **Geometric**: a subspace is coarser than the matrix — it discards the
   lengths of the basis vectors (the total counts per column) and keeps only
   directions. This cancels sequence-length and sequencing-depth differences.

If `rank(M) = k`, then `[M]` is a point of the Grassmann manifold `Gr(k, N)`.

---

## 5. Grassmann manifold

### 5.1 Definition

$$
\mathrm{Gr}(k,N)=\{\,V\subseteq\mathbb{R}^N:\ V \text{ is a } k\text{-dimensional linear subspace}\,\}.
$$

### 5.2 Geodesic distance

Let `A, B` be orthonormal bases of `V, W`, with principal angles
`θ₁ ≤ … ≤ θ_k`. The geodesic distance is

$$
d_{\mathrm{geo}}(V,W)=\|\Theta\|_F=\sqrt{\theta_1^2+\cdots+\theta_k^2}.
$$

It is intrinsic: independent of how the subspace is embedded in the ambient
space. This is what the paper means by a distance that differs from the
Euclidean distance between k-mer frequency vectors.

---

## 6. Principal angles

### 6.1 Definition (recursive)

For `dim P = p ≤ dim Q = q`,

$$
\cos\theta_1=\max_{\substack{x\in P,\ \|x\|=1\\ y\in Q,\ \|y\|=1}} x^{\mathsf T}y
=:x_1^{\mathsf T}y_1,
$$

$$
\cos\theta_k=\max_{\substack{x\in P,\ \|x\|=1,\ x^{\mathsf T}x_i=0\\ y\in Q,\ \|y\|=1,\ y^{\mathsf T}y_i=0}} x^{\mathsf T}y
=:x_k^{\mathsf T}y_k,\quad i<k .
$$

Geometrically: find the smallest angle and the corresponding pair of principal
vectors, then repeat on the orthogonal complements.

### 6.2 SVD computation (Proposition 1)

If `M_P ∈ R^{N×p}` and `M_Q ∈ R^{N×q}` have orthonormal columns and

$$
M_P^{\mathsf T}M_Q=U\Sigma V^{\mathsf T},\qquad
\Sigma=\begin{pmatrix}\tilde\Sigma&0\\0&0\end{pmatrix},\quad
\tilde\Sigma=\mathrm{diag}(\sigma_1,\dots,\sigma_p),
$$

then `θ_k = arccos(σ_k)`, with principal vectors `x_k = M_P u_k` and
`y_k = M_Q v_k`.

*Proof sketch.* For unit vectors `x = M_P u`, `y = M_Q v`,

$$
x^{\mathsf T}y=\sum_{i=1}^{p}\sigma_i\,(u^{\mathsf T}u_i)(v_i^{\mathsf T}v)
\le \sigma_1\sum_{i}|u^{\mathsf T}u_i|\,|v_i^{\mathsf T}v|\le\sigma_1,
$$

with equality at `u = u₁, v = v₁`; the recursive step is analogous. ∎

### 6.3 Sine form (Proposition 2)

$$
\bigl(I-M_QM_Q^{\mathsf T}\bigr)M_P \ \text{has singular values}\ \mu_k=\sqrt{1-\sigma_k^2}=\sin\theta_k .
$$

*Proof.* Let `A = (I − M_QM_Qᵀ)M_P`. Then

$$
A^{\mathsf T}A=I_p-(M_Q^{\mathsf T}M_P)^{\mathsf T}(M_Q^{\mathsf T}M_P)=I_p-U\Sigma\Sigma^{\mathsf T}U^{\mathsf T},
$$

so `UᵀAᵀAU = I_p − Σ̃²`, i.e. `μ_k = √(1−σ_k²)`. ∎

### 6.4 Numerical conditioning and the hybrid algorithm

`arccos(σ)` has condition number `1/√(1−σ²)` as `σ → 1`: when `θ ≈ 10⁻⁸`, `σ`
is within machine precision of 1 and `arccos` can lose all accuracy. `arcsin(μ)`
is well conditioned as `μ → 0` (condition number ≈ 1).

The implementation therefore uses a thresholded hybrid strategy with
`ε = cos(1e-5)` (small angle means `θ < 10⁻⁵` radians):

> **Algorithm 3 (principal angles)**
> Input: orthonormal bases `M_P ∈ R^{N×p}`, `M_Q ∈ R^{N×q}`, cosine threshold `ε`
> 1. Compute `M_PᵀM_Q = UΣVᵀ`, `Σ = diag(σ₁,…,σ_p)`
> 2. Compute the right principal vectors `V_cos = M_Q V`
> 3. If `σ₁ < ε`: return `θ_i = arccos(σ_i)` (all angles are "large")
> 4. Otherwise let `s = min{i : σ_i < ε}` (if none, `s = p+1`)
> 5. For `i ≥ s`: `θ_i = arccos(σ_i)`
> 6. Let `R_Q = V_cos[:, 1:s−1]`, compute `A = R_Q − M_P(M_PᵀR_Q)`
> 7. SVD of `A` gives `μ₁…μ_{s−1}`
> 8. For `i < s`: `θ_i = arcsin(μ_i)`
> 9. Return `θ` (ascending)

`fcgrsd.grassmann.principal_angles(method="hybrid")` implements this;
`method="acos"` is a pure-arccos variant for comparison.

---

## 7. Generalized Grassmann distance

The numerical ranks of different FCGRs usually differ (`k₁ ≠ k₂`), so the two
points live on Grassmannians of different dimensions and the geodesic distance
no longer applies. The paper uses the **generalized Grassmann distance** of
Ye & Lim (2016): it reduces to the geodesic distance when the dimensions match,
and otherwise equals the distance from one subspace to a particular Schubert
variety in the other Grassmannian.

Numerically it is still given by the principal angles:

$$
d(A_1,A_2)=\Bigl(\sum_{i=1}^{\min(k_1,k_2)}\theta_i^2\Bigr)^{1/2},
$$

where the number of angles is the **smaller** of the two dimensions.

> **Note (important).** Because only `min(k₁,k₂)` angles are used, this quantity
> is **not a metric** when the dimensions differ — the triangle inequality can
> fail. The smallest counterexample is in `R³`:
>
> $$A=\mathrm{span}(e_1)\subset B=\mathrm{span}(e_1,e_2)\supset C=\mathrm{span}(e_2).$$
>
> `A ⊆ B` and `C ⊆ B` give `d(A,B) = d(B,C) = 0`, but `A ⊥ C` gives
> `d(A,C) = π/2`, so `d(A,C) = 1.571 > d(A,B) + d(B,C) = 0`.
>
> In a random sample of 2,000 triples of subspaces in `R⁴⁰` (dimensions 1–19),
> about 6% violate the triangle inequality.
>
> This does not prevent NJ (which only needs a dissimilarity), but branch
> lengths should not be read as metric distances, and pairs with very different
> dimensions can be systematically compressed or stretched.

---

## 8. Truncation: why `p` is needed

`rank(M)` is well defined in theory, but in practice:

- floating-point SVD never returns exact zeros;
- sequencing errors introduce k-mers that do not exist in the genome;
- real FCGR spectra decay slowly, so the "true rank" is ambiguous.

The paper defines the numerical rank by the cumulative singular-value energy
`p ∈ (0,1)`:

$$
k=\#\Bigl\{i:\ \sum_{j\le i}d_j \le p\sum_{j}d_j\Bigr\}
\quad\text{or}\quad
k=\min\Bigl\{i:\ \sum_{j\le i}d_j \ge p\sum_{j}d_j\Bigr\}.
$$

The first is the literal wording of the paper (largest `k` with cumulative sum
not exceeding the target); the second is the smallest `k` that reaches the
target. They usually differ by one dimension. The package exposes both as
`rule="le"` and `rule="ge"` (default). The default reproduces the published
numbers; the difference in tree topology is usually small.

**Why truncation is principled.** By the Eckart–Young theorem, the truncated SVD
is the optimal low-rank approximation for a given rank, so keeping the first `k`
singular vectors is the best `k`-dimensional summary of the FCGR.

**Choosing `p`** (paper §1.2.4):

- small datasets: scan `0.01–0.99` in steps of `0.01` and take the value that
  minimizes the RF distance to a reference tree;
- large datasets: choose `p` on a random subset and reuse it.

`scan_truncation_p()` and `python -m fcgrsd scan-p` implement this protocol.

---

## 9. Full algorithm

```
Input: sequences S = {s_1,...,s_m}, resolution n, threshold p
Output: distance matrix D (m×m) and an NJ tree

for i = 1..m:
    M_i ← FCGR(s_i, n)                       # O(n·L_i)
    U_i, d_i, _ ← SVD(M_i)                   # O(8^n)
    k_i ← truncation rank from the energy rule
    A_i ← U_i[:, 1:k_i]                      # orthonormal basis, O(4^n·k_i)

for i < j:
    U, σ, V ← SVD(A_iᵀ A_j)                  # O(k_i k_j min(k_i,k_j))
    B ← (I - A_i A_iᵀ) A_j V                 # small-angle projection, O(4^n·k_i k_j)
    θ ← hybrid acos/asin principal angles
    D[i,j] ← ‖θ‖₂

T ← neighbor_joining(D)                      # O(m³)
return D, T
```

**Complexity**

| Stage | Complexity | Note |
|---|---|---|
| FCGR | `O(n·L)` | vectorized; can be reduced to `O(L)` with a rolling window |
| Truncated SVD | `O(m·8^n)` | full SVD of an `N=2^n` square matrix is `O(N³)` |
| Pairwise distance | `O(m²·(4^n·k² + k³))` | the projection `(I−A_iA_iᵀ)A_jV` is `O(N·k²)`; the angle SVD is `O(k³)`. For typical `n=6–8, k≈25`, the former dominates |
| NJ | `O(m³)` | standard cost of distance-based tree building |
| **Total** | `O(m·8^n + m²(4^n k² + k³) + m³)` | favorable for large `m` compared with MSA's `O(m²L²)` |

**Memory**: one FCGR is `O(4^n)` (`n=8` → 256×256); all truncated bases are
`O(m·4^n·k)`. The main pipeline keeps only the bases.

---

## 10. Parameter cheat sheet

| Parameter | Meaning | Recommended | Effect |
|---|---|---|---|
| `n` | FCGR resolution, matrix `2^n×2^n` | `ceil(log₄ L)`; HA≈6, SARS-CoV-2≈8 | too small loses detail; too large makes the matrix sparse and SVD slow |
| `p` | cumulative singular-value energy | scan on small data; 0.71 for HA, 0.75 for SARS-CoV-2 | too small underfits; too large admits noise dimensions |
| `alphabet` | base-to-corner convention (1–6) | 3 (paper default) | changes FCGR row/column structure; effect is limited |
| `rule` | truncation criterion | `ge` (default) or `le` | differ by one dimension; topology usually unchanged |

---

## 11. Properties and limitations

**Strengths**

- no alignment, no equal-length requirement, no normalization;
- only SVD and small matrix products, easy to parallelize;
- the distance is intrinsic to the manifold and insensitive to global scaling
  and length differences.

**Limitations**

1. **Resolution vs. sparsity.** When `4^n ≫ L`, the FCGR is extremely sparse and
   the column space is dominated by sampling noise. Keep `n ≈ log₄ L`.
2. **Distance saturation.** Unrelated sequences have nearly orthogonal subspaces,
   so distances approach the upper bound `√k·π/2`; highly divergent sequences
   are "flattened" and NJ loses resolution in deep branches.
3. **Not a metric across dimensions.** Taking only `min(k₁,k₂)` angles can
   violate the triangle inequality (see §7). NJ still works, but branch lengths
   are not metric distances.
4. **Parameter sensitivity.** The optimal `p` depends on dataset size and
   divergence; scan it rather than assuming a universal value.
