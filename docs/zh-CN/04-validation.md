# 验证报告

[English](../04-validation.md) | 简体中文

---

## 1. 验证策略

正确性以**独立参照**校验，而不是只与某一个实现比对：

- CGR 网格：`fractions.Fraction` 有理数精确递推；
- 主角度：`scipy.linalg.subspace_angles`；
- 截断维数与 Robinson-Foulds：暴力枚举；
- 邻接法拓扑：Biopython 的 NJ；
- MDS：欧氏构型重建。

目标是同时捕获数学错误与实现疏漏（下标偏移、顺序错误、dtype 静默提升）。

## 2. 测试套件

```bash
python -m pytest tests -q     # 113 passed
```

| 文件 | 用例数 | 覆盖内容 |
|---|---|---|
| `test_cgr.py` | 42 | 网格定位 vs 有理数精确 CGR（6 种字母表 × 5 种分辨率）；浮点/整数双路径；退化同聚物；坐标范围；异常输入 |
| `test_fcgr.py` | 16 | 计数守恒；单 k-mer 落格公式；归一化；非法碱基剔除与警告；大小写等价；`fcgr_total` |
| `test_grassmann.py` | 15 | 主角度 vs SciPy；小角度精度与升序；三角不等式；旋转不变性；截断单调性；参数校验 |
| `test_tree.py` | 15 | 可加矩阵下 NJ 精确还原；Newick 往返；已知 RF；引号名；二分划计数；NaN/重名拒绝；负枝长截断 |
| `test_pipeline.py` | 10 | 端到端结构；相同序列距离为 0；距离单调性；`ge` vs `le`；MDS 重建；重名/单序列拒绝 |
| `test_io.py` | 5 | 距离矩阵往返；表头与首列一致；方阵校验；空序列与重名 FASTA |
| `test_scripts.py` | 6 | 5 种登录号表格式的列识别 |
| `test_viz.py` | 2 | FCGR 热图朝向与 CGR 一致；全零枝长不产生退化坐标轴 |
| `test_cli.py` | 2 | `mds --dims > m` 的表头；`tree --clamp-negative` |

## 3. 数值校验

### CGR 网格

对 6 种字母表与 `n ∈ {1,2,3,6,8}`，把整数闭式解与"有理数精确递推 + `floor(2^n x)+1`"
逐元素比对，全部相同。浮点路径在随机序列与同聚物序列上与整数路径一致。

### 主角度

在随机子空间（含 4 个 `1e-7` 量级的小角度）上，与
`scipy.linalg.subspace_angles` 的最大差为 `9.7e-17`，且返回值升序。

### 距离

维数为 4 与 6 的两个子空间，`grassmann_distance` 与 SciPy 主角度的 ℓ₂ 范数
相差小于 `1e-12`。

### 截断

对 1000 组随机奇异值谱与随机 `p`，`truncation_rank` 与暴力扫描在 `ge`、`le`
两种规则下均一致。

### 邻接法

- 200 棵随机可加树：cophenetic 矩阵恢复到 `3.6e-15`，与真值树的 RF 为 0；
- 100 组随机非可加欧氏距离矩阵：拓扑与 Biopython 的 NJ 全部一致；
- 串行与 4 线程的距离矩阵逐位相同。

### Robinson-Foulds

对 100 组随机树对，基于二分划的 RF 与暴力对称差完全一致。

### MDS

对随机欧氏构型，经典 MDS 的距离重建误差为 `6.2e-15`。

## 4. 30 条 HA 数据上的复现

数据：`examples/data/tableS1_ha30.fasta`（30 条甲型流感 HA 基因，登录号见论文 Table S1）。

```bash
python -m fcgrsd run --fasta examples/data/tableS1_ha30.fasta \
    --k 6 --p 0.71 --out-prefix out/ha30
```

| 指标 | 结果 |
|---|---|
| 最优 `p`（0.01–0.99 扫描，RF 准则） | **0.71** |
| 截断维数 | 24–25（均值 ≈ 25） |
| 单系亚型 | **6 / 6** |
| 运行时间 | **0.039 s** |
| 与 NW+NJ 参考树的 RF | **4** |

6 个亚型（H1N1、H2N2、H3N2、H5N1、H7N3、H7N9）各自构成单系群，复现了论文
Figure 1B 的定性结论。论文报告的是相对 MAFFT + 最大似然参考树的 RF = 8；
本仓库的参考树用 Needleman–Wunsch + NJ 构建，绝对数值不可直接比较。

p 扫描结果在 `examples/results/tableS1_scan.csv`，曲线见
`examples/figures/fig4_scan_p.png`。

## 5. 模拟数据

`examples/data/demo_ha30.*` 包含 30 条 1700 bp 的模拟序列，沿已知随机树按
Jukes–Cantor 模型生成。以 `n=6, p=0.71` 运行，与真值树的 RF = 8（归一化 0.148），
与深度分化模拟枝上分辨率下降的预期一致。

```bash
python scripts/make_demo_data.py --out-dir examples/data --seed 20240909
python -m fcgrsd run --fasta examples/data/demo_ha30.fasta \
    --k 6 --p 0.71 --out-prefix out/demo
```

## 6. 性能与内存

30 条 HA 数据（`n=6`，`p=0.71`，单核）实测：

| 阶段 | 耗时 |
|---|---|
| FCGR + 截断 SVD | 0.015 s |
| 435 次成对距离 | 0.027 s |
| NJ | 0.0005 s |
| **合计** | **0.039 s** |

内存：

- `run`：同时只保留一个 FCGR 与一个完整 `U`，最终只存截断基（`O(m·4^n·k)`）；
- `scan-p`：保留所有序列的完整 `U`，`O(m·4^n)`（`n=8`、10348 条时约 5.4 GB），
  面向几十条序列。

## 7. 局限

以下是方法本身的性质，不是实现缺陷，记录在此以免过度解读：

1. **跨维数距离不是度量。** 只取 `min(k₁,k₂)` 个主角度时，三角不等式可能被违反。
   NJ 只要求相异性；枝长不是度量距离。
2. **`p` 的选择依赖数据。** 在同一批数据上以 RF 最小选 `p` 再报告 RF，会偏乐观。
3. **距离饱和。** 高度分化的序列子空间近似正交，距离趋近 `√k·π/2`，深层分支分辨率下降。
4. **非法碱基。** 剔除长 `N` 段会拼接两侧序列，可能产生伪 k-mer。
5. **没有置信值。** 流程输出点估计（一棵树），不提供 bootstrap 支持率。

## 8. 复现命令

```bash
# 安装
pip install -r requirements.txt

# 单元测试
python -m pytest tests -q

# 30 条 HA
python -m fcgrsd run --fasta examples/data/tableS1_ha30.fasta \
    --k 6 --p 0.71 --out-prefix out/ha30

# 对照 NW+NJ 参考树扫描 p
python -m fcgrsd scan-p --fasta examples/data/tableS1_ha30.fasta --k 6 \
    --reference examples/results/tableS1_scan_ref.nwk --out out/scan.csv

# 生成图（需 matplotlib）
python scripts/make_figures.py \
    --fasta examples/data/tableS1_ha30.fasta \
    --meta examples/data/tableS1_ha30_meta.csv \
    --k 6 --p 0.71 \
    --scan examples/results/tableS1_scan.csv \
    --out-dir examples/figures
```

> 参考树 `examples/results/tableS1_scan_ref.nwk` 是用于 p 扫描的
> Needleman–Wunsch + NJ 树。若更换输入数据，可用
> `scripts/benchmark_alignment.py` 重新生成。
