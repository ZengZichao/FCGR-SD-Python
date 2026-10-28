"""端到端流程与 MDS 的测试。"""

from __future__ import annotations

import numpy as np
import pytest

from fcgrsd.metrics import classical_mds
from fcgrsd.pipeline import run_fcgr_sd, suggest_n


def random_seq(n: int, seed: int = 0) -> str:
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list("ACGT"), n))


def mutate(seq: str, rate: float, seed: int) -> str:
    rng = np.random.default_rng(seed)
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


@pytest.fixture(scope="module")
def small_dataset():
    root = random_seq(1200, seed=1)
    names = [f"t{i}" for i in range(6)]
    seqs = [root] + [mutate(root, 0.02 * (i + 1), seed=10 + i) for i in range(5)]
    return names, seqs


def test_run_result_structure(small_dataset) -> None:
    names, seqs = small_dataset
    res = run_fcgr_sd(names, seqs, 6, 0.71)
    m = len(names)
    assert res.distances.shape == (m, m)
    assert np.allclose(np.diag(res.distances), 0.0)
    assert np.allclose(res.distances, res.distances.T)
    assert len(res.bases) == m
    assert res.ranks.shape == (m,)
    assert res.tree is not None and res.tree.n_leaves == m
    assert set(res.tree.leaf_names()) == set(names)
    for key in ("step1_2_fcgr_svd", "step3_distance", "step4_nj", "total"):
        assert key in res.timings


def test_identical_sequences_have_zero_distance(small_dataset) -> None:
    names, seqs = small_dataset
    res = run_fcgr_sd([names[0], names[0] + "_copy"], [seqs[0], seqs[0]], 5, 0.7)
    assert res.distances[0, 1] == pytest.approx(0.0, abs=1e-9)


def test_closer_sequences_have_smaller_distance(small_dataset) -> None:
    """与祖先序列的差异越小，Grassmann 距离应越小。"""
    names, seqs = small_dataset
    res = run_fcgr_sd(names, seqs, 6, 0.71)
    d = res.distances[0]
    assert d[1] < d[5]
    assert d[1] < d[3]


def test_ge_and_le_rules_are_close(small_dataset) -> None:
    """两种截断规则（差 1 维）应给出几乎相同的距离与完全相同的树拓扑。"""
    from fcgrsd.metrics import robinson_foulds

    names, seqs = small_dataset
    a = run_fcgr_sd(names, seqs, 6, 0.71, rule="ge")
    b = run_fcgr_sd(names, seqs, 6, 0.71, rule="le")
    assert np.max(np.abs(a.ranks - b.ranks)) <= 1
    upper = np.triu_indices(len(names), 1)
    corr = np.corrcoef(a.distances[upper], b.distances[upper])[0, 1]
    assert corr > 0.99
    assert robinson_foulds(a.tree, b.tree) == 0


def test_suggest_n() -> None:
    assert suggest_n(1700) == 6
    assert suggest_n(29759) == 8
    assert suggest_n(4) == 1


def test_duplicate_names_raise(small_dataset) -> None:
    names, seqs = small_dataset
    with pytest.raises(ValueError, match="唯一"):
        run_fcgr_sd([names[0], names[0]], seqs[:2], 5, 0.7)


def test_single_sequence_raises() -> None:
    with pytest.raises(ValueError, match="两条"):
        run_fcgr_sd(["A"], ["ACGTACGTAC"], 2, 0.7)


def test_mds_recovers_euclidean_configuration() -> None:
    rng = np.random.default_rng(3)
    X = rng.normal(size=(20, 3))
    D = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    coords, vals = classical_mds(D, 3)
    D2 = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1))
    assert np.allclose(D, D2, atol=1e-8)
    assert np.all(vals[:3] > 0)


def test_mds_handles_two_points() -> None:
    D = np.array([[0.0, 1.5], [1.5, 0.0]])
    coords, _ = classical_mds(D, 2)
    assert coords.shape == (2, 2)


def test_mds_clamps_components_to_point_count() -> None:
    """请求维数大于点数时应返回 m 列，而不是静默返回少于 n_components 的列。"""
    D = np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 1.5], [2.0, 1.5, 0.0]])
    coords, vals = classical_mds(D, 5)
    assert coords.shape == (3, 3)
    assert vals.shape == (3,)
