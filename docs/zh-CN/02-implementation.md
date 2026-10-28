# 实现说明

本文档说明包的组织方式与数值细节。数学推导见 [01 — 算法原理](01-algorithm.md)。

[English](../02-implementation.md) | 简体中文

---

## 1. 包结构

| 模块 | 职责 |
|---|---|
| `fcgrsd.cgr` | CGR 坐标、序列清洗、精确/浮点网格定位 |
| `fcgrsd.fcgr` | FCGR 构造、归一化、非法碱基警告 |
| `fcgrsd.grassmann` | 截断 SVD、主角度、广义距离 |
| `fcgrsd.distance` | 成对距离矩阵、可选线程并行、SVD 缓存 |
| `fcgrsd.tree` | 邻接法、Newick 读写、树结构查询 |
| `fcgrsd.metrics` | Robinson-Foulds 距离、经典 MDS |
| `fcgrsd.io` | FASTA 与距离矩阵读写 |
| `fcgrsd.pipeline` | 端到端流程与 p 扫描 |
| `fcgrsd.viz` | CGR、FCGR、树与 MDS 绘图 |
| `fcgrsd.cli` | `python -m fcgrsd` 入口 |

依赖关系无环：`cgr → fcgr → grassmann → distance → pipeline`，
`tree`/`metrics`/`io`/`viz` 由 pipeline 与 CLI 使用。

## 2. 数据流

```
FASTA ──read_fasta──▶ names, sequences
                         │
                         ▼
              fcgr_matrix(seq, n, alphabet)
                         │  （int64 计数矩阵，2^n × 2^n）
                         ▼
              subspace_basis(M, p, rule)      # 每条序列一次 SVD
                         │  （N × k 正交基）
                         ▼
        distance_matrix_from_bases(bases)     # m(m-1)/2 次小 SVD
                         │
                         ▼
              neighbor_joining(names, D) ──▶ Tree ──▶ Newick
```

`run_fcgr_sd` 同时只保留一个 FCGR 与一个完整 `U`，最终只存截断基，
因此峰值内存为 `O(4^n + m·4^n·k)`，而不是 `O(m·4^n)`。

## 3. 精确网格定位

由 CGR 递推可得格子编号的闭式解（推导见 [01](01-algorithm.md#23-格子编号的闭式解)）：

```
cell = 1 + Σ_{k=0}^{n-1} bit_{i-n+1+k} · 2^k
```

实现把每个碱基编码为两个比特向量 `u`、`v`，用 `n` 次向量化移位累加窗口和：

```python
windows = len(seq) - n + 1
rows = np.ones(windows, dtype=np.int64)
for k in range(n):
    rows += u[k : k + windows] * (1 << k)
```

复杂度为每条序列 `O(n·L)`，但循环体是单条向量操作，常数很小。结果是精确的：
不涉及浮点比较，同一个 k-mer 无论出现在序列的哪个位置，都映射到同一个格子。

`cgr_grid_indices_float` 用浮点递推实现同一映射（`ceil(2^n x)`，裁剪到 `[1, 2^n]`），
仅用于数值对照；所有入口默认 `mode="exact"`。

## 4. FCGR 构造

`fcgr_matrix` 把格子编号展平，用一次 `np.bincount` 完成计数：

```python
flat = (idx[:, 0] - 1) * size + (idx[:, 1] - 1)
counts = np.bincount(flat, minlength=size * size)
mat = counts.reshape(size, size)
```

**计数与频率**：论文称其为 frequency matrix。除以总数只是正数倍缩放，列空间不变，
距离完全相同。`normalize=True` 供需要概率的用户使用，主流程使用计数。

**非法碱基**：`normalize_sequence` 剔除 `A/C/G/T` 之外的字符并转大写。若一条序列中
非法字符超过 1%，会给出带具体占比的警告。剔除 `N` 是草图基因组的合理默认，
但长 `N` 段两侧会被直接拼接，可能产生基因组中不存在的 k-mer；必要时请按 `N` 切分。

**分辨率**：`suggest_n(L) = ceil(log₄ L)` 让格子数与序列长度同量级。它只是便捷函数，
`n` 可以显式指定。

## 5. 截断与 SVD 复用

截断维数由累积奇异值贡献率定义：

```python
cum = np.cumsum(s)
target = p * s.sum()
k = np.searchsorted(cum, target, side="left") + 1   # rule="ge"（默认）
# k = np.searchsorted(cum, target, side="right")     # rule="le"
```

`ge` 取达到目标的最小 `k`；`le` 取不超过目标的最大 `k`。两者至多相差 1 维。
默认 `ge` 以复现论文数值；`le` 对应论文文字的字面表述。返回的秩会夹取到 `[1, N]`。

**SVD 缓存**：p 扫描需要对 99 个 `p` 复用同一次 SVD。`build_svd_cache` 为每条序列
保存 `(U, s)`，`basis_from_svd` 按不同 `p` 切片。这是用内存换时间：缓存保留完整 `U`，
即 `O(m·4^n)`（`n=8`、10348 条时约 5.4 GB），只适合论文规模的数据集（30–44 条）。
主流程 `run` **不**使用该缓存。

## 6. 主角度

`principal_angles` 计算 `A₁ᵀA₂ = UΣVᵀ`，并以 `ε = cos(1e-5)` 分段：

- `σ ≥ ε` 的小角度用投影矩阵 `B = (I − A₁A₁ᵀ)A₂V` 的奇异值 `μ` 做 `arcsin`；
- 其余角度用 `arccos(σ)`。

`σ` 与 `μ` 在反三角函数前都裁剪到 `[0, 1]`。舍入可能给出 `σ = 1 + 9e-16`；
不裁剪时 `arccos` 返回复数，整个距离矩阵会被提升为复数类型。

小角度块在 `arcsin` 前按升序排序，因此返回向量满足 `θ₁ ≤ θ₂ ≤ …`。距离是范数，
对置换不变，但顺序属于公开契约，并有回归测试覆盖。

## 7. 距离矩阵

`distance_matrix_from_bases` 只填上三角并镜像：

```python
for i in range(m - 1):
    for j in range(i + 1, m):
        d = grassmann_distance(bases[i], bases[j])
        D[i, j] = D[j, i] = d
```

每个无序对只算一次，因此结果严格对称。`n_jobs > 1` 时按行分配到
`ThreadPoolExecutor`；LAPACK 在 SVD 期间释放 GIL，因此能获得真实并行，
且基矩阵无需序列化。各线程写入 `D` 的不同单元。

## 8. 邻接法与 Newick

实现的是经典 Saitou–Nei 算法：

```
while 剩余结点多于两个:
    r_i = Σ_k D_ik
    Q_ij = (n-2) D_ij - r_i - r_j
    合并使 Q 最小的结点对
    d(i,u) = D_ij/2 + (r_i - r_j)/(2(n-2))
    d(j,u) = D_ij - d(i,u)
    D_uk = (D_ik + D_jk - D_ij)/2
```

距离矩阵原地更新；失效的行/列置零并用掩码排除，避免复制子矩阵。最后两个结点各分
剩余距离的一半，形成二叉根。并列最优时取展平行优先 `argmin` 的最小者，结果确定可复现。

距离矩阵不可加时 NJ 可能给出负枝长。默认保留；`clamp_negative=True`
（CLI 的 `--clamp-negative`）会截断为 0。

Newick 输出会给含空白或 `()[],:;'"` 的标签加引号，并去掉枝长末尾的零。
`read_newick` 支持引号标签、`[...]` 注释与枝长。

## 9. Robinson-Foulds 与 MDS

`robinson_foulds` 默认比较无根树的二分划，`rooted=True` 时比较有根树的 clade。
平凡二分划（单叶、或除一叶外的全部）被忽略。归一化值为
`|S_a Δ S_b| / (|S_a| + |S_b|)`，落在 `[0, 1]`。叶名重复的树会被拒绝：
二分划集合会塌缩。

`classical_mds` 是 Torgerson 的 PCoA：

```python
B = -0.5 * J @ (D**2) @ J
vals, vecs = np.linalg.eigh(B)
coords = vecs[:, :k] * np.sqrt(np.clip(vals[:k], 0, None))
```

负特征值截断为 0，这是非欧氏相异性做 PCoA 的标准处理。返回的维数为
`min(n_components, m)`。

## 10. 错误处理

| 输入 | 行为 |
|---|---|
| 序列长度小于 `n` | `ValueError` |
| 空序列或全为非法字符 | `ValueError` |
| `n` 不是正整数 | `ValueError` |
| `p` 不在 `(0, 1)` | `ValueError` |
| `alphabet` 不在 `1..6` | `ValueError` |
| 序列重名 | `ValueError`（FASTA 读取、pipeline、NJ） |
| 距离矩阵含 `NaN`/`Inf` | `ValueError` |
| 距离矩阵表头与首列标签不一致 | `ValueError` |
| `run` 少于两条序列 | `ValueError` |
| `mds --dims` 大于点数 | 夹取到 `m` |

目标是：对会产出"看起来合理但没有意义"的树，直接报错。

## 11. 内存与性能

- `run_fcgr_sd`：同时只保留一个 FCGR 与一个完整 `U`，最终只存截断基
  （`O(m·4^n·k)`）。
- `scan_truncation_p`：保留所有序列的完整 `U`（`O(m·4^n)`），以避免 99 次重复 SVD。
- 成对距离：输出 `O(m²)`；每对一次 `k³` 的 SVD 加一次 `N·k²` 的投影。
- 并行：`--n-jobs` 并行成对循环；使用线程而非进程，基矩阵无需序列化。

在 30 条 HA 数据上实测（`n=6`，`p=0.71`）：端到端约 0.04 s，平均截断维数约 25。
