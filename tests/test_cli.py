"""CLI 子命令的回归测试。"""

from __future__ import annotations

import csv
from types import SimpleNamespace

import numpy as np

from fcgrsd.cli import cmd_mds, cmd_tree
from fcgrsd.io import write_distance_matrix
from fcgrsd.tree import read_newick


def _lengths(node):
    out = []
    if node.length is not None:
        out.append(node.length)
    for c in node.children:
        out.extend(_lengths(c))
    return out


def test_cmd_mds_header_matches_columns(workdir) -> None:
    """--dims 大于点数时，表头列数必须与实际坐标列数一致。"""
    D = np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 1.5], [2.0, 1.5, 0.0]])
    dist = workdir / "d.csv"
    out = workdir / "mds.csv"
    write_distance_matrix(str(dist), ["A", "B", "C"], D)
    cmd_mds(SimpleNamespace(dist=str(dist), dims=5, out=str(out), plot=None))
    rows = list(csv.reader(open(out, encoding="utf-8")))
    assert len(rows[0]) == len(rows[1]) == 4  # name + min(5, 3) 维


def test_cmd_tree_clamp_negative(workdir) -> None:
    D = np.array(
        [
            [0.0, 2.5782, 1.8943, 1.0337],
            [2.5782, 0.0, 1.7420, 0.6057],
            [1.8943, 1.7420, 0.0, 0.3070],
            [1.0337, 0.6057, 0.3070, 0.0],
        ]
    )
    dist = workdir / "d.csv"
    out = workdir / "t.nwk"
    write_distance_matrix(str(dist), ["L0", "L1", "L2", "L3"], D)
    cmd_tree(SimpleNamespace(dist=str(dist), out=str(out), clamp_negative=True))
    tree = read_newick(out.read_text(encoding="utf-8"))
    assert min(_lengths(tree.root)) >= 0
