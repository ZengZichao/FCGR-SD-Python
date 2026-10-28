# 示例数据

[English](README.md) | 简体中文

## `data/`

| 文件 | 说明 |
|---|---|
| `tableS1_ha30.fasta` | 从 NCBI 下载的 30 条甲型流感 HA 基因（登录号见论文 Table S1） |
| `tableS1_ha30_meta.csv` | `name,subtype` 元数据 |
| `demo_ha30.fasta` | 30 条 1700 bp 模拟序列，带已知真值树 |
| `demo_ha30_meta.csv` | 模拟分组标签 |
| `demo_ha30_true.nwk` | 模拟数据的真值 Newick 树 |

## `results/`

| 文件 | 说明 |
|---|---|
| `tableS1_scan.csv` | p 扫描结果（99 个 p 值 + 2 条 FCGR 基线），用于 `fig4_scan_p.png` |
| `tableS1_scan_ref.nwk` | 扫描用的 Needleman–Wunsch + NJ 参考树 |

## `figures/`

由 `scripts/make_figures.py` 生成：

- `fig1_cgr_fcgr.png` —— CGR 走笔与两种分辨率的 FCGR；
- `fig2_tree.png` —— 30 条 HA 的 NJ 树，按亚型着色；
- `fig3_mds.png` —— 距离矩阵的经典 MDS；
- `fig4_scan_p.png` —— RF 与平均截断维数随 `p` 的变化。
