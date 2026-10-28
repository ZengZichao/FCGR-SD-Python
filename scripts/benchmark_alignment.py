"""Reproduce the RF-based evaluation protocol of the paper (Table 1).

The paper uses MAFFT + maximum-likelihood trees as the reference. To avoid an
external MSA dependency, this script uses Biopython's Needleman-Wunsch global
alignment to build pairwise p-distances and then the same NJ builder as the
main method. The absolute RF values will differ from the paper (different
reference tree), but the "scan p, minimize RF" protocol of §1.2.4 is the same,
so the script is useful for parameter selection and relative comparison.

用法::

    python scripts/benchmark_alignment.py --fasta examples/data/tableS1_ha30.fasta \
        --k 6 --meta examples/data/tableS1_ha30_meta.csv --out examples/results/tableS1_scan.csv
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Dict, List, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fcgrsd.distance import build_bases, build_svd_cache, distance_matrix_from_bases  # noqa: E402
from fcgrsd.fcgr import fcgr_matrix  # noqa: E402
from fcgrsd.io import read_fasta  # noqa: E402
from fcgrsd.metrics import robinson_foulds  # noqa: E402
from fcgrsd.pipeline import scan_truncation_p  # noqa: E402
from fcgrsd.tree import neighbor_joining  # noqa: E402


def alignment_distances(seqs: List[str]) -> np.ndarray:
    """全局比对（NW）得到的 p-距离矩阵：1 - 比对一致率。"""
    from Bio import Align

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -1.0

    m = len(seqs)
    D = np.zeros((m, m))
    for i in range(m - 1):
        for j in range(i + 1, m):
            aln = aligner.align(seqs[i], seqs[j])[0]
            counts = aln.counts()
            ident = counts.identities / max(counts.aligned, 1)
            D[i, j] = D[j, i] = 1.0 - ident
    return D


def fcgr_baseline_distances(seqs: List[str], k: int, metric: str = "euclidean") -> np.ndarray:
    """FCGR 基线（Hatje & Kollmar 的 FCGR 方法）：直接比较频率矩阵。"""
    mats = [fcgr_matrix(s, k, 3, normalize=True).ravel() for s in seqs]
    X = np.vstack(mats)
    if metric == "euclidean":
        # 逐行计算，避免 m×m×4^k 的三维临时数组（k=8、m=30 时约 470 MB）
        D = np.zeros((X.shape[0], X.shape[0]))
        for i in range(X.shape[0]):
            D[i] = np.sqrt(((X - X[i]) ** 2).sum(axis=1))
        return D
    if metric == "pearson":
        Xc = X - X.mean(axis=1, keepdims=True)
        norm = np.linalg.norm(Xc, axis=1, keepdims=True)
        norm[norm == 0] = 1.0  # 常数向量（理论上不会出现）避免除零
        Xc = Xc / norm
        return 1.0 - Xc @ Xc.T
    raise ValueError(metric)


def clade_purity(tree, names: List[str], meta: Dict[str, str]) -> Tuple[int, int]:
    """统计有多少个分组在树上表现为单系群。"""
    clusters = tree.clusters()
    all_names = set(names)
    groups: Dict[str, set] = {}
    for n in names:
        groups.setdefault(meta.get(n, ""), set()).add(n)
    ok = 0
    for g, s in groups.items():
        if not g:
            continue
        if any(s <= c for c in clusters) or any(s <= (all_names - c) for c in clusters):
            ok += 1
    return ok, len([g for g in groups if g])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--meta", default=None)
    ap.add_argument("--out", default="scan_p.csv")
    ap.add_argument("--p-grid", default="0.01:0.99:0.01")
    args = ap.parse_args()

    names, seqs = read_fasta(args.fasta)
    print(f"序列数 {len(names)}，长度中位数 {int(np.median([len(s) for s in seqs]))}")

    meta: Dict[str, str] = {}
    if args.meta and os.path.exists(args.meta):
        with open(args.meta, encoding="utf-8") as fh:
            meta = {r["name"]: r["group"] for r in csv.DictReader(fh)}

    print("计算全局比对距离矩阵（参考树）...")
    D_aln = alignment_distances(seqs)
    ref_tree = neighbor_joining(names, D_aln)
    with open(os.path.splitext(args.out)[0] + "_ref.nwk", "wt", encoding="utf-8") as fh:
        fh.write(ref_tree.to_newick() + "\n")

    lo, hi, step = (float(x) for x in args.p_grid.split(":"))
    ps = np.round(np.arange(lo, hi + 1e-9, step), 4)

    print("扫描截断参数 p ...")
    rows = scan_truncation_p(names, seqs, args.k, ps=ps, reference_tree=ref_tree, verbose=True)

    # FCGR 基线
    for metric in ("euclidean", "pearson"):
        D = fcgr_baseline_distances(seqs, args.k, metric)
        t = neighbor_joining(names, D)
        rows.append(
            {
                "p": float("nan"),
                "rank_mean": float("nan"),
                "rank_min": 0,
                "rank_max": 0,
                "rf": robinson_foulds(t, ref_tree),
                "rf_normalized": robinson_foulds(t, ref_tree, normalize=True),
                "note": f"FCGR-{metric}",
            }
        )

    best = min((r for r in rows if r["p"] == r["p"]), key=lambda r: (r["rf"], r["p"]))
    print(
        f"\nRF 最小的 FCGR-SD 参数：p = {best['p']}，RF = {best['rf']:.0f}，"
        f"平均截断维数 {best['rank_mean']:.1f}"
    )

    # 在最优 p 与论文推荐 p 上检查单系性
    cache = build_svd_cache(seqs, args.k)
    for p in sorted({best["p"], 0.71}):
        bases = build_bases(seqs, args.k, p, cache=cache)
        tree = neighbor_joining(names, distance_matrix_from_bases(bases))
        if meta:
            ok, tot = clade_purity(tree, names, meta)
            print(f"  p={p}: 单系群 {ok}/{tot}，RF={robinson_foulds(tree, ref_tree):.0f}")
        with open(f"{os.path.splitext(args.out)[0]}_p{p}.nwk", "wt", encoding="utf-8") as fh:
            fh.write(tree.to_newick() + "\n")

    with open(args.out, "wt", encoding="utf-8", newline="") as fh:
        keys = ["p", "rank_mean", "rank_min", "rank_max", "rf", "rf_normalized", "note"]
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"扫描结果：{args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
