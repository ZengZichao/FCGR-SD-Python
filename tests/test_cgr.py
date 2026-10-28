"""CGR 坐标与网格定位的单元测试。"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

from fcgrsd.cgr import (
    ALPHABETS,
    cgr_coordinates,
    cgr_grid_indices,
    cgr_grid_indices_float,
    encode_bits,
    normalize_sequence,
)


def random_seq(n: int, seed: int = 0) -> str:
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list("ACGT"), n))


def exact_cells(seq: str, n: int, kap: int) -> np.ndarray:
    """用有理数精确递推 CGR，再按 floor(2^n x)+1 定位（独立参照实现）。"""
    mapping = ALPHABETS[kap]
    x = Fraction(1, 2)
    y = Fraction(1, 2)
    pts = []
    for ch in seq:
        u, v = mapping[ch]
        x = (x + u) / 2
        y = (y + v) / 2
        pts.append((x, y))
    size = 2**n
    return np.array([[int(px * size) + 1, int(py * size) + 1] for px, py in pts[n - 1 :]])


@pytest.mark.parametrize("kap", [1, 2, 3, 4, 5, 6])
@pytest.mark.parametrize("n", [1, 2, 3, 6, 8])
def test_grid_indices_match_exact_rational(kap: int, n: int) -> None:
    """整数快速路径必须与有理数精确递推逐元素一致（全部 6 种字母表）。"""
    seq = random_seq(400, seed=kap * 10 + n)
    assert np.array_equal(cgr_grid_indices(seq, n, kap), exact_cells(seq, n, kap))


@pytest.mark.parametrize("kap", [1, 3, 6])
def test_float_path_matches_exact(kap: int) -> None:
    """浮点路径与整数路径在常规序列上应完全一致。"""
    seq = random_seq(3000, seed=7)
    n = 8
    assert np.array_equal(cgr_grid_indices(seq, n, kap), cgr_grid_indices_float(seq, n, kap))


def test_float_path_handles_degenerate_runs() -> None:
    """长同聚物（浮点退化场景）下两条路径仍应一致。"""
    seq = "ACGT" + "G" * 500 + "ACGT" + "A" * 500 + "ACGT" + "C" * 200
    n = 6
    assert np.array_equal(cgr_grid_indices(seq, n, 3), cgr_grid_indices_float(seq, n, 3))


def test_coordinates_in_unit_square() -> None:
    seq = random_seq(500, seed=3)
    coords = cgr_coordinates(seq, 3)
    assert coords.shape == (len(seq), 2)
    assert np.all(coords >= 0.0) and np.all(coords <= 1.0)


def test_first_point_is_half_corner() -> None:
    """首个 CGR 点应为 (1/2 + 顶点)/2，即位于对应象限中心。"""
    assert np.allclose(cgr_coordinates("A", 3)[0], [0.25, 0.25])
    assert np.allclose(cgr_coordinates("C", 3)[0], [0.75, 0.75])
    assert np.allclose(cgr_coordinates("G", 3)[0], [0.25, 0.75])
    assert np.allclose(cgr_coordinates("T", 3)[0], [0.75, 0.25])


def test_exact_coordinates_match_fraction() -> None:
    seq = random_seq(60, seed=11)
    approx = cgr_coordinates(seq, 3)
    exact = cgr_coordinates(seq, 3, exact=True)
    # 精确模式无累积误差；两种模式在短序列上应当一致到 1e-15
    assert np.allclose(approx, exact, atol=1e-15)


def test_alphabet_permutation_changes_geometry() -> None:
    """不同 ka 必须给出不同的 CGR 几何（否则说明查表未生效）。"""
    seq = random_seq(200, seed=5)
    a = cgr_coordinates(seq, 1)
    b = cgr_coordinates(seq, 3)
    assert not np.allclose(a, b)


def test_normalize_sequence_filters_and_uppercases() -> None:
    assert normalize_sequence("acgtNNac-gt") == "ACGTACGT"
    with pytest.raises(ValueError):
        normalize_sequence("ACGTN", on_unknown="error")


def test_encode_bits_shape_and_values() -> None:
    u, v = encode_bits("ACGT", 3)
    assert u.shape == v.shape == (4,)
    assert set(np.unique(u)).issubset({0, 1})


def test_short_sequence_raises() -> None:
    with pytest.raises(ValueError):
        cgr_grid_indices("ACG", 6, 3)


def test_invalid_alphabet_raises() -> None:
    with pytest.raises(ValueError):
        cgr_grid_indices("ACGTACGT", 2, 7)
