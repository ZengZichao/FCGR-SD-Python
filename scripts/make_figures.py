"""生成文档与示例用的图表（CGR/FCGR、进化树、MDS、p 扫描曲线）。

用法::

    python scripts/make_figures.py --fasta examples/data/tableS1_ha30.fasta \
        --meta examples/data/tableS1_ha30_meta.csv --k 6 --p 0.71 \
        --scan examples/results/tableS1_scan.csv --out-dir examples/figures
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Dict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from fcgrsd.fcgr import fcgr_matrix  # noqa: E402
from fcgrsd.io import read_fasta  # noqa: E402
from fcgrsd.metrics import classical_mds  # noqa: E402
from fcgrsd.pipeline import run_fcgr_sd  # noqa: E402
from fcgrsd.viz import draw_tree, plot_cgr, plot_fcgr, plot_mds  # noqa: E402

PALETTE = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
]


def group_colors(names, meta: Dict[str, str]) -> Dict[str, str]:
    groups = sorted({meta.get(n, "") for n in names if meta.get(n, "")})
    cmap = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(groups)}
    return {n: cmap.get(meta.get(n, ""), "#333333") for n in names}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", default="examples/data/tableS1_ha30.fasta")
    ap.add_argument("--meta", default="examples/data/tableS1_ha30_meta.csv")
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--p", type=float, default=0.71)
    ap.add_argument("--scan", default="examples/results/tableS1_scan.csv")
    ap.add_argument("--out-dir", default="examples/figures")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    names, seqs = read_fasta(args.fasta)
    meta: Dict[str, str] = {}
    if os.path.exists(args.meta):
        with open(args.meta, encoding="utf-8") as fh:
            meta = {r["name"]: r["group"] for r in csv.DictReader(fh)}
    colors = group_colors(names, meta) if meta else {}

    # --- 图 1：CGR 走笔 + 两种分辨率的 FCGR ---
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.8), dpi=170)
    plot_cgr(seqs[0], axes[0], walk=300)
    for ax, k in ((axes[1], 4), (axes[2], args.k)):
        im = plot_fcgr(fcgr_matrix(seqs[0], k, 3), ax)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, label="log(1+count)")
    fig.suptitle(f"{names[0]}  (L={len(seqs[0])})", fontsize=9)
    fig.tight_layout()
    f1 = os.path.join(args.out_dir, "fig1_cgr_fcgr.png")
    fig.savefig(f1)
    plt.close(fig)
    print("已写出", f1)

    # --- 图 2：进化树 ---
    res = run_fcgr_sd(names, seqs, args.k, args.p)
    fig, ax = plt.subplots(figsize=(6.6, 5.6), dpi=180)
    draw_tree(res.tree, ax, colors=colors, label_size=6.2)
    ax.set_title(f"FCGR-SD NJ tree  (n={args.k}, p={args.p})", fontsize=9)
    fig.tight_layout()
    f2 = os.path.join(args.out_dir, "fig2_tree.png")
    fig.savefig(f2)
    plt.close(fig)
    print("已写出", f2)

    # --- 图 3：MDS ---
    coords, vals = classical_mds(res.distances, 3)
    fig = plt.figure(figsize=(6.0, 5.2), dpi=180)
    ax = fig.add_subplot(111, projection="3d")
    plot_mds(coords, names, ax, colors=colors, dims=3, s=30)
    ax.set_title("Classical MDS of Grassmann distances", fontsize=9)
    fig.tight_layout()
    f3 = os.path.join(args.out_dir, "fig3_mds.png")
    fig.savefig(f3)
    plt.close(fig)
    print("已写出", f3)

    # --- 图 4：p 扫描曲线 ---
    if os.path.exists(args.scan):
        ps, rfs, ranks = [], [], []
        with open(args.scan, encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                if not row["p"] or row["p"] != row["p"] or row.get("note"):
                    continue
                try:
                    ps.append(float(row["p"]))
                    rfs.append(float(row["rf"]))
                    ranks.append(float(row["rank_mean"]))
                except ValueError:
                    continue
        fig, ax1 = plt.subplots(figsize=(5.6, 3.6), dpi=180)
        ax1.plot(ps, rfs, color="#c53030", lw=1.4, label="RF distance")
        ax1.set_xlabel("truncation parameter p", fontsize=8)
        ax1.set_ylabel("RF distance", color="#c53030", fontsize=8)
        ax1.tick_params(labelsize=7)
        ax2 = ax1.twinx()
        ax2.plot(ps, ranks, color="#2b6cb0", lw=1.0, ls="--", label="mean rank")
        ax2.set_ylabel("mean truncation rank", color="#2b6cb0", fontsize=8)
        ax2.tick_params(labelsize=7)
        best = int(np.argmin(rfs))
        ax1.scatter([ps[best]], [rfs[best]], color="#c53030", zorder=5, s=22)
        ax1.annotate(
            f"p*={ps[best]:.2f}",
            (ps[best], rfs[best]),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=7,
            color="#c53030",
        )
        ax1.set_title("RF vs truncation parameter", fontsize=9)
        fig.tight_layout()
        f4 = os.path.join(args.out_dir, "fig4_scan_p.png")
        fig.savefig(f4)
        plt.close(fig)
        print("已写出", f4)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
