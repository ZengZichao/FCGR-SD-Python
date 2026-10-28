"""Generate a reproducible simulated dataset with a known true tree.

The repository ships the 30 HA sequences from NCBI, but a simulated dataset
with a known topology is useful for quantifying reconstruction accuracy (RF).
The script builds a random binary tree, simulates Jukes-Cantor substitutions
along its branches, and writes FASTA, the true Newick tree, and group metadata.

Usage::

    python scripts/make_demo_data.py --out-dir examples/data --n-groups 6 --length 1700
"""

from __future__ import annotations

import argparse
import os
import random
from dataclasses import dataclass, field
from typing import List, Optional

BASES = "ACGT"


@dataclass
class SimNode:
    name: Optional[str] = None
    length: float = 0.0
    children: List["SimNode"] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children


def random_group(names: List[str], rng: random.Random, lo: float, hi: float) -> SimNode:
    """把一组名字随机组装成二叉子树（枝长落在 [lo, hi]）。"""
    pool = [SimNode(name=n, length=rng.uniform(lo, hi)) for n in names]
    while len(pool) > 1:
        i = rng.randrange(len(pool))
        a = pool.pop(i)
        j = rng.randrange(len(pool))
        b = pool.pop(j)
        pool.append(SimNode(length=rng.uniform(lo, hi), children=[a, b]))
    return pool[0]


def mutate(seq: str, t: float, rng: random.Random) -> str:
    """Jukes-Cantor 模拟：每个位点以 p = 3/4 (1 - e^{-4t/3}) 的概率突变。"""
    import math

    p = 0.75 * (1.0 - math.exp(-4.0 * t / 3.0))
    out = list(seq)
    for i, ch in enumerate(out):
        if rng.random() < p:
            others = [b for b in BASES if b != ch]
            out[i] = rng.choice(others)
    return "".join(out)


def simulate(node: SimNode, seq: str, rng: random.Random, acc: dict) -> None:
    if node.is_leaf:
        acc[node.name] = seq
        return
    for child in node.children:
        simulate(child, mutate(seq, child.length, rng), rng, acc)


def to_newick(node: SimNode, digits: int = 6) -> str:
    if node.is_leaf:
        return f"{node.name}:{node.length:.{digits}f}"
    inner = ",".join(to_newick(c, digits) for c in node.children)
    if node.length > 0:
        return f"({inner}):{node.length:.{digits}f}"
    return f"({inner})"


def main() -> None:
    ap = argparse.ArgumentParser(description="生成 FCGR-SD 演示数据")
    ap.add_argument("--out-dir", default="examples/data")
    ap.add_argument("--n-groups", type=int, default=6)
    ap.add_argument("--per-group", type=int, default=5)
    ap.add_argument("--length", type=int, default=1700)
    ap.add_argument("--seed", type=int, default=20240909)
    ap.add_argument("--name", default="demo_ha30")
    ap.add_argument(
        "--within",
        type=float,
        nargs=2,
        default=(0.002, 0.012),
        metavar=("LO", "HI"),
        help="类群内部枝长区间（对应同亚型 ~1-3%% 差异）",
    )
    ap.add_argument(
        "--between",
        type=float,
        nargs=2,
        default=(0.06, 0.12),
        metavar=("LO", "HI"),
        help="类群之间（根部到类群祖先）的枝长区间",
    )
    args = ap.parse_args()

    rng = random.Random(args.seed)
    root_seq = "".join(rng.choice(BASES) for _ in range(args.length))

    groups = []
    all_names = []
    for g in range(args.n_groups):
        names = [f"G{g + 1}_{i + 1:02d}" for i in range(args.per_group)]
        all_names.extend(names)
        groups.append((f"group{g + 1}", names))

    subtrees = [
        random_group(names, rng, args.within[0], args.within[1]) for _, names in groups
    ]
    for st in subtrees:
        st.length = rng.uniform(args.between[0], args.between[1])
    # 把 6 个类群随机组装成二叉根，避免真值树出现多分支（便于计算 RF）
    root = subtrees[0]
    for st in subtrees[1:]:
        root = SimNode(length=rng.uniform(*args.between), children=[root, st])
    root.length = 0.0

    acc: dict = {}
    simulate(root, root_seq, rng, acc)

    os.makedirs(args.out_dir, exist_ok=True)
    fa = os.path.join(args.out_dir, f"{args.name}.fasta")
    with open(fa, "wt", encoding="utf-8") as fh:
        for name in all_names:
            fh.write(f">{name}\n")
            s = acc[name]
            for i in range(0, len(s), 70):
                fh.write(s[i : i + 70] + "\n")

    nwk = os.path.join(args.out_dir, f"{args.name}_true.nwk")
    with open(nwk, "wt", encoding="utf-8") as fh:
        fh.write(to_newick(root) + ";\n")

    meta = os.path.join(args.out_dir, f"{args.name}_meta.csv")
    with open(meta, "wt", encoding="utf-8") as fh:
        fh.write("name,group\n")
        for gname, names in groups:
            for n in names:
                fh.write(f"{n},{gname}\n")

    print(f"已生成 {len(all_names)} 条序列（长度 {args.length}）")
    print(f"  FASTA  : {fa}")
    print(f"  真值树 : {nwk}")
    print(f"  元数据 : {meta}")


if __name__ == "__main__":
    main()
