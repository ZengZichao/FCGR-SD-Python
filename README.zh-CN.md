# FCGR-SD

**基于 Grassmann 流形的无比对基因组比较** —— 论文

> Li X., Zhou T., Feng X., Yau S.-T., Yau S. S.-T. (2024).
> *Exploring geometry of genome space via Grassmann manifolds.*
> **The Innovation** 5(5): 100677. <https://doi.org/10.1016/j.xinn.2024.100677>

的 Python 实现。

[English](README.md) | 简体中文

---

## 方法

```
DNA 序列 ──CGR──▶ 单位正方形上的点集 ──网格计数──▶ FCGR 矩阵 M ∈ R^{N×N}
                                                        │
                                                        │ 取列空间
                                                        ▼
                                          Grassmann 流形 Gr(k, N) 上的点 [A]
                                                        │
                                                        │ 广义 Grassmann 距离
                                                        ▼
                                            d([A₁],[A₂]) ──NJ──▶ 系统发育树
```

1. **CGR**：把每个碱基映射到单位正方形的一个顶点，当前点向该顶点移动一半。
2. **FCGR**：统计落入 `2^n × 2^n` 网格各格子的 CGR 点数。
3. **截断 SVD**：保留前若干个左奇异向量，即列空间的一组标准正交基 `A`。
4. **主角度**：由 `A₁ᵀA₂` 的 SVD 得到两个列空间之间的主角度。
5. **广义 Grassmann 距离**：主角度向量的 ℓ₂ 范数（支持维数不同的子空间）。
6. **邻接法（NJ）**：由距离矩阵构建系统发育树。

无需多序列比对、不要求序列等长、不直接比较 k-mer 频率向量，只依赖 SVD 与少量矩阵乘法。

## 特点

- **整数精确网格定位**：每个 k-mer 的 CGR 格子由闭式解直接算出，不存在浮点格线平局。
- **向量化 FCGR + 截断 SVD**：每条序列只做一次 SVD，只保留截断基。
- **混合 `arccos`/`arcsin` 主角度**：对几乎相同的子空间仍保持精度。
- **完整工具链**：NJ、Newick 读写、Robinson-Foulds 距离、经典 MDS、p 扫描、命令行。
- **113 项单元测试**，使用独立参照（`scipy`、Biopython、暴力枚举、有理数精确递推）。

## 安装

```bash
git clone https://github.com/ZengZichao/FCGR-SD-Python.git
cd FCGR-SD-Python
python -m venv .venv && source .venv/bin/activate
pip install -e .            # 核心依赖：numpy
# 或：pip install -r requirements.txt   # 含测试/绘图/基准的可选依赖
```

Python ≥ 3.9。核心依赖 `numpy`；可选 `scipy`（测试）、`matplotlib`（绘图）、
`biopython`（数据下载与比对基准）、`pytest`。

## 快速开始

```bash
python -m fcgrsd run \
    --fasta examples/data/tableS1_ha30.fasta \
    --k 6 --p 0.71 \
    --out-prefix out/ha30
```

输出：

- `out/ha30.dist.csv` —— 带表头与行标签的距离矩阵；
- `out/ha30.nwk` —— Newick 格式的邻接树；
- `out/ha30.json` —— 运行摘要（参数、截断维数、距离统计、耗时）。

## 命令行

| 命令 | 用途 |
|---|---|
| `fcgrsd run` | FASTA → 距离矩阵 → NJ 树 |
| `fcgrsd scan-p` | 扫描截断阈值 `p`（可对照参考树） |
| `fcgrsd tree` | 由已有距离矩阵构建 NJ 树 |
| `fcgrsd mds` | 对距离矩阵做经典 MDS 降维 |

常用选项：`--k/--n`、`--p`、`--alphabet 1..6`、`--rule ge|le`、
`--mode exact|float`、`--n-jobs`、`--clamp-negative`、`--quiet`。
详见 `python -m fcgrsd <command> --help`。

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

底层构件：

```python
from fcgrsd.fcgr import fcgr_matrix
from fcgrsd.grassmann import subspace_basis, grassmann_distance

M1 = fcgr_matrix(seqs[0], n=6)
M2 = fcgr_matrix(seqs[1], n=6)
A1 = subspace_basis(M1, p=0.71)
A2 = subspace_basis(M2, p=0.71)
print(grassmann_distance(A1, A2))
```

## 复现 30 条 HA 的例子

仓库附带论文 Figure 1B 使用的 30 条甲型流感 HA 基因（登录号见论文 Table S1，来源于 NCBI）。

```bash
python -m fcgrsd run --fasta examples/data/tableS1_ha30.fasta \
    --k 6 --p 0.71 --out-prefix out/ha30
```

该数据集上的预期结果：

| 指标 | 数值 |
|---|---|
| 最优截断阈值 | `p = 0.71`（网格扫描，RF 准则） |
| 平均截断维数 | ≈ 25（范围 24–25） |
| 单系亚型 | 6 / 6 |
| 运行时间（30 条） | ≈ 0.04 s |
| 与 NW+NJ 参考树的 RF | 4 |

论文报告的是相对 MAFFT + 最大似然参考树的 RF = 8；本仓库的参考树用
Needleman–Wunsch + NJ 构建，绝对 RF 不可直接比较。定性结论——6 个亚型各自成单系群——已复现。

## 目录结构

```
FCGR-SD-Python/
├── fcgrsd/                 Python 包
│   ├── cgr.py              CGR 坐标与精确网格定位
│   ├── fcgr.py             FCGR 构造
│   ├── grassmann.py        截断 SVD、主角度、广义距离
│   ├── distance.py         成对距离矩阵
│   ├── tree.py             邻接法与 Newick 读写
│   ├── metrics.py          Robinson-Foulds 距离与经典 MDS
│   ├── io.py               FASTA 与距离矩阵读写
│   ├── pipeline.py         端到端流程
│   ├── viz.py              CGR/FCGR/树/MDS 绘图
│   └── cli.py              命令行入口
├── tests/                  113 项单元测试
├── scripts/                数据下载、模拟、绘图、基准
├── examples/
│   ├── data/               30 条 HA、元数据、模拟数据
│   ├── results/            预计算的 p 扫描结果
│   └── figures/            文档插图
├── docs/                   英文文档
│   └── zh-CN/              中文文档
├── pyproject.toml
├── requirements.txt
├── CITATION.cff
└── LICENSE
```

## 文档

| 文档 | 内容 |
|---|---|
| [01 — 算法原理](docs/zh-CN/01-algorithm.md) | CGR/FCGR 定义、Grassmann 几何、主角度、距离、复杂度 |
| [02 — 实现说明](docs/zh-CN/02-implementation.md) | 包设计、精确网格推导、数值细节、错误处理 |
| [03 — 使用手册](docs/zh-CN/03-usage.md) | 安装、数据准备、CLI 与 API、参数调优、FAQ |
| [04 — 验证报告](docs/zh-CN/04-validation.md) | 独立校验、单元测试、复现结果、性能与局限 |

英文版本：[`docs/`](docs/)。

## 验证

核心实现以独立参照校验，而不是只与某一个实现比对：

| 模块 | 参照 | 结果 |
|---|---|---|
| 网格定位 | 有理数精确 CGR（`Fraction`） | 逐元素一致 |
| 主角度 | `scipy.linalg.subspace_angles` | 最大差 9.7e-17 |
| 距离 | SciPy 主角度的 ℓ₂ 范数 | 差 < 1e-12 |
| 截断维数 | 暴力枚举 | 1000 组谱 0 不一致 |
| NJ | Biopython `DistanceTreeConstructor.nj` | 100 组随机矩阵拓扑 0 不一致 |
| RF | 暴力枚举二分划 | 100 棵树 0 不一致 |
| MDS | 欧氏构型重建 | 最大误差 6.2e-15 |

运行测试：`python -m pytest tests -q`。

## 局限

- **跨维数距离不是度量。** 只取 `min(k₁,k₂)` 个主角度时，三角不等式可能被违反
  （可能出现 `d(A,C) > d(A,B) + d(B,C)`）。NJ 只要求相异性，建树仍可进行，
  但枝长不应解释为度量距离。
- **`p` 依赖数据。** 论文通过最小化与参考树的 RF 来选 `p`；若在同一批数据上
  既调参又评估，结果会偏乐观。验证性分析建议使用留出子集或固定 `p`。
- **非法碱基会被剔除。** 这对草图基因组是合理默认，但剔除长 `N` 段会把两侧序列
  直接拼接，可能产生基因组中并不存在的 k-mer；必要时请按 `N` 切分。
- **可能出现负枝长。** 距离矩阵不可加时 NJ 会给出负枝长；默认保留，
  需要非负时用 `--clamp-negative`。
- **`scan-p` 会缓存完整的左奇异矩阵**，以避免 99 次重复 SVD。内存为 `O(m·4^n)`
  （`n=8`、10348 条时约 5.4 GB），只适合论文规模的数据集（30–44 条）。
  主流程 `run` 是流式的，只保留截断基。

## 数据

- `examples/data/tableS1_ha30.fasta` —— 从 NCBI 下载的 30 条甲型流感 HA 基因
  （公共数据；登录号见论文 Table S1）。
- `examples/data/demo_ha30.*` —— 由 `scripts/make_demo_data.py` 生成的模拟数据，
  带有已知真值树。
- 仓库**不重新分发** GISAID 数据。`scripts/fetch_datasets.py` 可按你自己的登录号表
  从 NCBI 下载；GISAID 登录号会被记为缺失。

## 引用

如果使用本软件，请引用论文：

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

## 许可

MIT，见 [LICENSE](LICENSE)。
