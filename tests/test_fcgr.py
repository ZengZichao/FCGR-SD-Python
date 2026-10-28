"""FCGR 矩阵的单元测试。"""

from __future__ import annotations

import numpy as np
import pytest

from fcgrsd.cgr import ALPHABETS
from fcgrsd.fcgr import fcgr_matrix


def random_seq(n: int, seed: int = 0) -> str:
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list("ACGT"), n))


@pytest.mark.parametrize("n", [1, 2, 4, 6])
def test_total_count(n: int) -> None:
    """FCGR 总计数必须等于参与定位的 CGR 点数 L-n+1。"""
    seq = random_seq(1000, seed=n)
    M = fcgr_matrix(seq, n, 3)
    assert M.shape == (2**n, 2**n)
    assert int(M.sum()) == len(seq) - n + 1


def test_single_kmer_cell_formula() -> None:
    """单条 n-mer 的落格位置应满足 1 + sum_k bit_k 2^k（首碱基权重 2^0）。"""
    n = 4
    for kmer in ("AAAA", "CCCC", "GGGG", "TTTT", "ACGT", "TGCA", "AGCT", "GATC"):
        M = fcgr_matrix(kmer, n, 3)
        assert M.sum() == 1
        row, col = np.argwhere(M == 1)[0]
        mp = ALPHABETS[3]
        exp_r = 1 + sum(mp[ch][0] * (1 << k) for k, ch in enumerate(kmer))
        exp_c = 1 + sum(mp[ch][1] * (1 << k) for k, ch in enumerate(kmer))
        assert (row + 1, col + 1) == (exp_r, exp_c)


def test_float_mode_matches_exact_mode() -> None:
    seq = random_seq(1500, seed=21)
    assert np.array_equal(fcgr_matrix(seq, 6, 3, mode="exact"), fcgr_matrix(seq, 6, 3, mode="float"))


def test_normalize_sums_to_one() -> None:
    seq = random_seq(500, seed=2)
    M = fcgr_matrix(seq, 5, 3, normalize=True)
    assert M.dtype == np.float64
    assert pytest.approx(M.sum(), rel=1e-12) == 1.0


def test_unknown_bases_are_skipped() -> None:
    base = random_seq(300, seed=4)
    noisy = base[:100] + "N" + base[100:] + "-"
    assert np.array_equal(fcgr_matrix(base, 5, 3), fcgr_matrix(noisy, 5, 3))


def test_lowercase_equivalence() -> None:
    seq = random_seq(300, seed=6)
    assert np.array_equal(fcgr_matrix(seq, 5, 3), fcgr_matrix(seq.lower(), 5, 3))


def test_short_sequence_raises() -> None:
    with pytest.raises(ValueError):
        fcgr_matrix("ACGT", 6, 3)


def test_ambiguous_fraction() -> None:
    from fcgrsd.cgr import ambiguous_fraction

    assert ambiguous_fraction("ACGTACGT") == 0.0
    assert ambiguous_fraction("ACGTNNNN") == pytest.approx(0.5)
    assert ambiguous_fraction("") == 0.0


def test_warns_when_many_ambiguous_bases() -> None:
    """非 ACGT 占比 >1% 时必须显式提示，避免静默丢数据。"""
    import fcgrsd.fcgr as F

    with pytest.warns(UserWarning, match="非 ACGT"):
        F.fcgr_matrix("ACGT" * 50 + "N" * 20, 4, 3)


def test_warns_for_each_distinct_fraction() -> None:
    """不同严重程度的序列应各自提示，而不是第一条之后全部静默。"""
    import warnings

    import fcgrsd.fcgr as F

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        F.fcgr_matrix("ACGT" * 50 + "N" * 50, 4, 3)   # 50/250 = 20%
        F.fcgr_matrix("ACGT" * 5 + "N" * 80, 4, 3)    # 80/100 = 80%
    assert len(caught) == 2
    assert "80.0%" in str(caught[1].message)


def test_no_warning_on_clean_sequence() -> None:
    import warnings

    import fcgrsd.fcgr as F

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        F.fcgr_matrix("ACGT" * 100, 4, 3)


def test_fcgr_total_ignores_ambiguous_bases() -> None:
    from fcgrsd.fcgr import fcgr_total

    assert fcgr_total("ACGTACGT", 2) == 7
    assert fcgr_total("ACGTNNNN", 2) == 3          # 只数 4 个 ACGT
    assert fcgr_total("ACGT-acgt", 2) == 7
    assert fcgr_total("NNNN", 2) == 0


def test_non_integer_n_raises() -> None:
    with pytest.raises(ValueError):
        fcgr_matrix("ACGTACGT", 2.0, 3)
    with pytest.raises(ValueError):
        fcgr_matrix("ACGTACGT", True, 3)
