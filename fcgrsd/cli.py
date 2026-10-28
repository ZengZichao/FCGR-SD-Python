"""命令行入口：``python -m fcgrsd ...``"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import List, Optional

import numpy as np

from . import __version__
from .grassmann import DEFAULT_P
from .io import read_distance_matrix, read_fasta, write_distance_matrix
from .metrics import classical_mds
from .pipeline import run_fcgr_sd, scan_truncation_p, suggest_n
from .tree import neighbor_joining, read_newick


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--k", "--n", dest="k", type=int, default=None, help="FCGR 分辨率参数 n（矩阵 2^n × 2^n）")
    parser.add_argument("--p", type=float, default=DEFAULT_P, help=f"奇异值累积贡献率截断阈值（默认 {DEFAULT_P}）")
    parser.add_argument("--alphabet", type=int, default=3, choices=range(1, 7), help="碱基-顶点对应编号 1~6（默认 3）")
    parser.add_argument("--rule", choices=["ge", "le"], default="ge", help="截断维数判定规则：ge（默认）或 le")
    parser.add_argument("--mode", choices=["exact", "float"], default="exact", help="FCGR 网格定位方式（默认 exact）")
    parser.add_argument("--n-jobs", type=int, default=1, help="成对距离并行线程数")
    parser.add_argument("--quiet", action="store_true", help="关闭进度输出")


def cmd_run(args: argparse.Namespace) -> int:
    names, seqs = read_fasta(args.fasta)
    k = args.k if args.k is not None else suggest_n(int(np.median([len(s) for s in seqs])))
    if args.out_prefix:
        prefix = args.out_prefix
        os.makedirs(os.path.dirname(os.path.abspath(prefix)) or ".", exist_ok=True)
    else:
        prefix = os.path.splitext(os.path.basename(args.fasta))[0]

    res = run_fcgr_sd(
        names,
        seqs,
        k,
        args.p,
        alphabet=args.alphabet,
        rule=args.rule,
        mode=args.mode,
        n_jobs=args.n_jobs,
        build_tree=True,
        clamp_negative=args.clamp_negative,
        verbose=not args.quiet,
    )

    write_distance_matrix(f"{prefix}.dist.csv", names, res.distances)
    with open(f"{prefix}.nwk", "wt", encoding="utf-8") as fh:
        fh.write(res.tree.to_newick() + "\n")
    with open(f"{prefix}.json", "wt", encoding="utf-8") as fh:
        json.dump(res.summary(), fh, ensure_ascii=False, indent=2)
    if not args.quiet:
        print(f"已写出：{prefix}.dist.csv / {prefix}.nwk / {prefix}.json")
    return 0


def cmd_scan_p(args: argparse.Namespace) -> int:
    names, seqs = read_fasta(args.fasta)
    k = args.k if args.k is not None else suggest_n(int(np.median([len(s) for s in seqs])))
    ref = None
    if args.reference:
        with open(args.reference, encoding="utf-8") as fh:
            ref = read_newick(fh.read())
    rows = scan_truncation_p(
        names,
        seqs,
        k,
        reference_tree=ref,
        alphabet=args.alphabet,
        rule=args.rule,
        mode=args.mode,
        n_jobs=args.n_jobs,
        verbose=not args.quiet,
    )
    out = args.out or "scan_p.csv"
    keys = list(rows[0].keys())
    with open(out, "wt", encoding="utf-8", newline="") as fh:
        import csv

        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    best = None
    if ref is not None:
        # RF 并列时取较小的 p（截断维数更低、模型更简约）
        best = min(rows, key=lambda r: (r["rf"], r["p"]))
        print(f"RF 最小的 p = {best['p']}（RF={best['rf']}, 归一化={best['rf_normalized']:.4f}）")
    print(f"扫描结果已写出：{out}")
    return 0


def cmd_tree(args: argparse.Namespace) -> int:
    names, D = read_distance_matrix(args.dist)
    tree = neighbor_joining(names, D, clamp_negative=args.clamp_negative)
    out = args.out or "tree.nwk"
    with open(out, "wt", encoding="utf-8") as fh:
        fh.write(tree.to_newick() + "\n")
    print(f"已写出系统发育树：{out}")
    return 0


def cmd_mds(args: argparse.Namespace) -> int:
    if args.dims < 1:
        print("--dims 必须为正整数", file=sys.stderr)
        return 2
    names, D = read_distance_matrix(args.dist)
    coords, vals = classical_mds(D, args.dims)
    out = args.out or "mds.csv"
    with open(out, "wt", encoding="utf-8", newline="") as fh:
        import csv

        w = csv.writer(fh)
        # 用 coords.shape[1]（= min(dims, m)）而非 args.dims，避免点数少于维数时列数错位
        w.writerow(["name"] + [f"dim{i + 1}" for i in range(coords.shape[1])])
        for nm, row in zip(names, coords):
            w.writerow([nm] + [f"{v:.6f}" for v in row])
    print(f"已写出 MDS 坐标：{out}（前 3 个特征值：{np.round(vals[:3], 4)}）")
    if args.plot:
        _plot_mds(coords, names, vals, args.plot)
        print(f"已写出散点图：{args.plot}")
    return 0


def _plot_mds(coords, names, vals, path: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6.4, 5.2), dpi=160)
    if coords.shape[1] >= 3:
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=28, c="#2b6cb0")
        for i, nm in enumerate(names):
            ax.text(coords[i, 0], coords[i, 1], coords[i, 2], nm, fontsize=6)
        ax.set_xlabel("MDS1")
        ax.set_ylabel("MDS2")
        ax.set_zlabel("MDS3")
    elif coords.shape[1] == 2:
        ax = fig.add_subplot(111)
        ax.scatter(coords[:, 0], coords[:, 1], s=28, c="#2b6cb0")
        for i, nm in enumerate(names):
            ax.text(coords[i, 0], coords[i, 1], nm, fontsize=6)
        ax.set_xlabel("MDS1")
        ax.set_ylabel("MDS2")
    else:
        # 只保留 1 维时退化为数轴上的点
        ax = fig.add_subplot(111)
        ax.scatter(coords[:, 0], np.zeros(len(coords)), s=28, c="#2b6cb0")
        for i, nm in enumerate(names):
            ax.text(coords[i, 0], 0.0, nm, fontsize=6)
        ax.set_xlabel("MDS1")
        ax.set_yticks([])
    ax.set_title("FCGR-SD classical MDS")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fcgrsd",
        description="FCGR-SD：基于 Grassmann 流形的无比对基因组比较（The Innovation 2024）",
    )
    p.add_argument("--version", action="version", version=f"fcgrsd {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="完整流程：FASTA → 距离矩阵 → NJ 树")
    run.add_argument("--fasta", required=True)
    run.add_argument("--out-prefix", default=None)
    _add_common(run)
    run.add_argument("--clamp-negative", action="store_true", help="把 NJ 的负枝长截断为 0")
    run.set_defaults(func=cmd_run)

    scan = sub.add_parser("scan-p", help="扫描截断参数 p（可选：以参考树 RF 为准则）")
    scan.add_argument("--fasta", required=True)
    scan.add_argument("--reference", default=None, help="参考树 Newick 文件")
    scan.add_argument("--out", default=None)
    _add_common(scan)
    scan.set_defaults(func=cmd_scan_p)

    tree = sub.add_parser("tree", help="由距离矩阵构建 NJ 树")
    tree.add_argument("--dist", required=True)
    tree.add_argument("--out", default=None)
    tree.add_argument("--clamp-negative", action="store_true", help="把负枝长截断为 0")
    tree.set_defaults(func=cmd_tree)

    mds = sub.add_parser("mds", help="对距离矩阵做经典 MDS 降维")
    mds.add_argument("--dist", required=True)
    mds.add_argument("--dims", type=int, default=3)
    mds.add_argument("--out", default=None)
    mds.add_argument("--plot", default=None)
    mds.set_defaults(func=cmd_mds)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
