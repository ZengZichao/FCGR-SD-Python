# 使用手册

[English](../03-usage.md) | 简体中文

---

## 1. 安装

```bash
git clone https://github.com/ZengZichao/FCGR-SD-Python.git
cd FCGR-SD-Python
python -m venv .venv && source .venv/bin/activate
pip install -e .
```

要求 Python ≥ 3.9、`numpy ≥ 1.22`。可选依赖：

| 用途 | 包 | 使用场景 |
|---|---|---|
| 测试 | `pytest`、`scipy` | `pytest tests` |
| 绘图 | `matplotlib` | `scripts/make_figures.py`、`fcgrsd mds --plot` |
| 基准 | `scipy`、`biopython` | `scripts/benchmark_alignment.py`、`scripts/fetch_datasets.py` |

一次性安装全部：`pip install -r requirements.txt`。

检查安装：

```bash
python -m pytest tests -q      # 期望：113 passed
python -m fcgrsd --help
```

## 2. 数据准备

### FASTA

每条序列一条记录，名字唯一，使用标准 IUPAC 字符。CGR 之前会剔除 `N`、`-` 等
非法字符；单条序列中占比超过 1% 时给出警告。

```
>KM409247.1
ATGAAGGCAATACTAGTAGTTCTGCTATACACATTTGCAACCACAAATGCAGACACATTATGTATAGGTT...
```

序列名取标题行第一个空白前的 token，因此 `>seq1 description` 的名字是 `seq1`。

### 登录号表

`scripts/fetch_datasets.py` 接受带表头的 CSV。登录号列按表头中的 `accession`、
`gisaid`、`header`、`id` 依次识别；`Header` 形如 `"LC710779.1 Avian coronavirus"`
时取首 token `LC710779.1`。可选的分组列按 `gene`、`type`、`variant`、`subgenus`、
`lineage`、`pangolin`、`clade` 识别。

```bash
python scripts/fetch_datasets.py \
    --table my_accessions.csv \
    --out data/sequences.fasta \
    --batch 50
```

GISAID 登录号（`EPI_ISL_*`）通常无法从 NCBI 获取，会被写入 `*_missing.txt`。

## 3. 命令行

### `run` —— 从 FASTA 到树

```bash
python -m fcgrsd run \
    --fasta examples/data/tableS1_ha30.fasta \
    --k 6 --p 0.71 \
    --out-prefix out/ha30
```

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

| 选项 | 默认 | 含义 |
|---|---|---|
| `--fasta` | 必填 | 输入 FASTA |
| `--k` / `--n` | `ceil(log₄ L)` | FCGR 分辨率，矩阵 `2^k × 2^k` |
| `--p` | 0.75 | 奇异值累积贡献率阈值 |
| `--alphabet` | 3 | 碱基-顶点约定，1–6 |
| `--rule` | `ge` | 截断判据：`ge` 或 `le` |
| `--mode` | `exact` | 网格定位：`exact` 或 `float` |
| `--n-jobs` | 1 | 成对距离的并行线程数 |
| `--clamp-negative` | 关 | 把 NJ 负枝长截断为 0 |
| `--out-prefix` | FASTA 文件名 | 输出前缀 |
| `--quiet` | 关 | 关闭进度输出 |

输出：

- `.dist.csv` —— 表头 + 行标签 + 10 位小数距离；
- `.nwk` —— Newick 树；
- `.json` —— 参数、截断维数统计、距离统计、耗时。

### `scan-p` —— 选择截断阈值

```bash
python -m fcgrsd scan-p \
    --fasta examples/data/tableS1_ha30.fasta \
    --k 6 \
    --reference examples/results/tableS1_scan_ref.nwk \
    --out out/scan.csv
```

在 `p ∈ [0.01, 0.99]` 上以 0.01 为步长扫描，每个 `p` 输出一行：
`p, rank_mean, rank_min, rank_max, rf, rf_normalized`。给出 `--reference` 时，
打印 RF 最小的 `p`（并列时取较小的 `p`）。

> `scan-p` 为每条序列缓存完整的 `U`，以避免 99 次重复 SVD。内存为 `O(m·4^n)`，
> 面向几十条序列，不适合上万条。

### `tree` —— 由距离矩阵建树

```bash
python -m fcgrsd tree --dist out/ha30.dist.csv --out out/ha30.nwk
python -m fcgrsd tree --dist out/ha30.dist.csv --out out/ha30.nwk --clamp-negative
```

适合距离矩阵由其它工具算出的场景。

### `mds` —— 经典多维标度

```bash
python -m fcgrsd mds --dist out/ha30.dist.csv --dims 3 \
    --out out/mds.csv --plot out/mds.png
```

`--dims` 大于点数时会夹取到 `m`，因此 CSV 表头始终与实际坐标列数一致。

## 4. Python API

### 端到端

```python
from fcgrsd.io import read_fasta
from fcgrsd.pipeline import run_fcgr_sd

names, seqs = read_fasta("examples/data/tableS1_ha30.fasta")
res = run_fcgr_sd(names, seqs, n=6, p=0.71)

res.distances          # (30, 30) 对称矩阵
res.ranks              # 每条序列的截断维数
res.tree               # Tree 对象
res.timings            # 各阶段耗时
res.summary()          # 可 JSON 序列化的摘要
```

`run_fcgr_sd` 支持 `rule="ge"|"le"`、`mode="exact"|"float"`、
`alphabet=1..6`、`n_jobs`、`build_tree`、`clamp_negative`。

### 底层构件

```python
from fcgrsd.fcgr import fcgr_matrix
from fcgrsd.grassmann import subspace_basis, principal_angles, grassmann_distance

M1 = fcgr_matrix(seqs[0], n=6)
M2 = fcgr_matrix(seqs[1], n=6)

A1 = subspace_basis(M1, p=0.71)          # (64, k1) 正交基
A2 = subspace_basis(M2, p=0.71)

theta = principal_angles(A1, A2)         # 升序主角度
d = grassmann_distance(A1, A2)           # ||theta||_2
```

### 树与度量

```python
from fcgrsd.tree import neighbor_joining, read_newick
from fcgrsd.metrics import robinson_foulds, classical_mds

tree = neighbor_joining(names, res.distances)
reference = read_newick(open("reference.nwk").read())
print(robinson_foulds(tree, reference))
print(robinson_foulds(tree, reference, normalize=True))

coords, eigvals = classical_mds(res.distances, 3)
```

### 扫描 p

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

## 5. 参数调优

### `n` —— 分辨率

让格子数 `4^n` 与序列长度同量级：

| 序列长度 | 建议 `n` | 矩阵 |
|---|---|---|
| < 500 | 4 | 16 × 16 |
| 1–3 kb（HA、spike） | 6 | 64 × 64 |
| 10–40 kb（SARS-CoV-2） | 8 | 256 × 256 |
| > 100 kb（细菌基因组） | 9–10 | 512²–1024² |

过小会丢失 k-mer 细节；过大会使矩阵稀疏、SVD 变慢。`n > 11` 会触发内存警告。

### `p` —— 截断

`p` 控制保留多少奇异值谱。论文推荐值为 30 条 HA 的 0.71、44 条 SARS-CoV-2 的 0.75。
新数据集应扫描 `p`（对照参考树或留出子集）。若把 RF 当作性能估计，不要在同一批数据上
既调参又评估。

### `alphabet`

`alphabet=3`（AGCT）是论文约定与默认值。其它约定会旋转/镜像 CGR 正方形，
可用于敏感性检查，但会改变 FCGR 的行列结构。

### `rule`

`ge`（默认）取累积贡献率达到 `p·total` 的最小 `k`；`le` 取不超过它的最大 `k`。
两者至多相差 1 维。只有需要严格按论文文字复现时才用 `le`。

### `mode`

`exact` 是默认值，生产环境应使用它。`float` 走浮点 CGR 递推，仅用于数值对照。

## 6. FAQ

**树里出现负枝长。**
距离矩阵不可加时 NJ 会给出负枝长，属正常现象。下游工具要求非负时用
`--clamp-negative`。

**序列含大量 `N`。**
默认剔除，占比超过 1% 时警告。若 `N` 成簇出现，剔除会把两侧拼接，可能产生伪 k-mer；
建议按 `N` 切分后分别计数。

**内存不够。**
主流程 `run` 是流式的。`scan-p` 内存不足时，减少序列数或降低 `n`；
p 扫描缓存为 `O(m·4^n)`。

**能否用于蛋白质序列？**
不能直接使用。方法依赖 DNA 的四顶点几何；蛋白质需要 20 顶点构造与不同的网格编码。

**RF 与论文不一致。**
论文以 MAFFT + 最大似然树为参考；本包的基准脚本用 Needleman–Wunsch + NJ。
绝对 RF 不可直接比较，但 p 选择流程与定性结论一致。

**为什么 `run` 至少需要两条序列？**
单条序列没有成对距离，也没有有意义的树。只需要 1×1 矩阵时，直接用
`distance_matrix`。

**能否使用自定义距离矩阵？**
可以：`fcgrsd tree --dist matrix.csv`、`fcgrsd mds --dist matrix.csv`。
CSV 需带表头、行标签，且为对称方阵。
