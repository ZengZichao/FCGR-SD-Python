"""NJ 建树、Newick 读写与 RF 距离的单元测试。"""

from __future__ import annotations

import numpy as np
import pytest

from fcgrsd.metrics import robinson_foulds
from fcgrsd.tree import neighbor_joining, read_newick


def test_nj_recovers_additive_tree_exactly() -> None:
    """可加距离矩阵下 NJ 应精确还原树度量（cophenetic 距离）。"""
    t = read_newick("((A:1,B:1):1,(C:1,D:1):1);")
    names, D = t.patristic_distances()
    rebuilt = neighbor_joining(names, D)
    names2, D2 = rebuilt.patristic_distances()
    order = [names2.index(n) for n in names]
    assert np.allclose(D, D2[np.ix_(order, order)], atol=1e-12)


def test_nj_topology_on_known_case() -> None:
    t = read_newick("((A:1,B:1):1,(C:1,D:1):1);")
    names, D = t.patristic_distances()
    rebuilt = neighbor_joining(names, D)
    assert robinson_foulds(rebuilt, t) == 0


def test_nj_two_taxa() -> None:
    tree = neighbor_joining(["A", "B"], np.array([[0.0, 2.0], [2.0, 0.0]]))
    names, D = tree.patristic_distances()
    assert names == ["A", "B"]
    assert D[0, 1] == pytest.approx(2.0)


def test_nj_single_taxon() -> None:
    tree = neighbor_joining(["A"], np.zeros((1, 1)))
    assert tree.n_leaves == 1


def test_nj_rejects_nan() -> None:
    with pytest.raises(ValueError, match="NaN"):
        neighbor_joining(["A", "B"], np.array([[0.0, np.nan], [np.nan, 0.0]]))


def test_nj_rejects_duplicate_names() -> None:
    with pytest.raises(ValueError, match="唯一"):
        neighbor_joining(["A", "A"], np.array([[0.0, 1.0], [1.0, 0.0]]))


def _all_lengths(node):
    out = []
    if node.length is not None:
        out.append(node.length)
    for c in node.children:
        out.extend(_all_lengths(c))
    return out


def test_nj_clamp_negative() -> None:
    """非可加距离会产生负枝长；clamp_negative=True 应把它们截断为 0。"""
    D = np.array(
        [
            [0.0, 2.5782, 1.8943, 1.0337],
            [2.5782, 0.0, 1.7420, 0.6057],
            [1.8943, 1.7420, 0.0, 0.3070],
            [1.0337, 0.6057, 0.3070, 0.0],
        ]
    )
    names = ["L0", "L1", "L2", "L3"]
    raw = neighbor_joining(names, D)
    assert min(_all_lengths(raw.root)) < 0
    clamped = neighbor_joining(names, D, clamp_negative=True)
    assert min(_all_lengths(clamped.root)) >= 0


def test_newick_roundtrip() -> None:
    text = "((A:0.1,B:0.2)90:0.3,(C:0.4,D:0.5):0.6);"
    tree = read_newick(text)
    assert tree.n_leaves == 4
    again = read_newick(tree.to_newick())
    assert robinson_foulds(tree, again) == 0
    assert again.leaf_names() == ["A", "B", "C", "D"]


def test_newick_quoted_names() -> None:
    tree = read_newick("('a b':1,'c,d':2);")
    assert set(tree.leaf_names()) == {"a b", "c,d"}
    out = tree.to_newick()
    assert "'a b'" in out


def test_rf_identical_trees() -> None:
    t1 = read_newick("((A:1,B:1):1,(C:1,D:1):1);")
    t2 = read_newick("((C:2,D:2):2,(A:2,B:2):2);")
    assert robinson_foulds(t1, t2) == 0


def test_rf_different_topology() -> None:
    t1 = read_newick("((A:1,B:1):1,(C:1,D:1):1);")
    t2 = read_newick("((A:1,C:1):1,(B:1,D:1):1);")
    assert robinson_foulds(t1, t2) == 2


def test_rf_normalized_range() -> None:
    t1 = read_newick("((A:1,B:1):1,(C:1,D:1):1);")
    t2 = read_newick("((A:1,C:1):1,(B:1,D:1):1);")
    val = robinson_foulds(t1, t2, normalize=True)
    assert 0.0 < val <= 1.0


def test_rf_requires_same_leaves() -> None:
    t1 = read_newick("((A:1,B:1):1,C:1);")
    t2 = read_newick("((A:1,D:1):1,C:1);")
    with pytest.raises(ValueError):
        robinson_foulds(t1, t2)


def test_patristic_distances_of_star_tree() -> None:
    t = read_newick("(A:1,B:1,C:1);")
    names, D = t.patristic_distances()
    assert names == ["A", "B", "C"]
    assert np.allclose(D, np.array([[0, 2, 2], [2, 0, 2], [2, 2, 0]]))


def test_splits_count_for_binary_tree() -> None:
    """m 个叶子的二叉无根树应有 m-3 个非平凡二分划。"""
    t = read_newick("(((A:1,B:1):1,C:1):1,(D:1,E:1):1);")
    assert len(t.splits()) == t.n_leaves - 3
