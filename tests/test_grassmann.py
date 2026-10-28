"""Grassmann 距离核心（截断 SVD、主角度、广义距离）的单元测试。"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import subspace_angles

from fcgrsd.grassmann import (
    DEFAULT_EPS,
    grassmann_distance,
    principal_angles,
    subspace_basis,
    truncation_rank,
)


def random_orthonormal(n: int, k: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return np.linalg.qr(rng.normal(size=(n, k)))[0]


def test_principal_angles_match_scipy() -> None:
    """主角度（升序）应与 scipy.linalg.subspace_angles（降序）一致。"""
    A1 = random_orthonormal(64, 5, seed=1)
    A2 = random_orthonormal(64, 7, seed=2)
    mine = principal_angles(A1, A2)
    ref = subspace_angles(A1, A2)
    assert np.allclose(mine, np.sort(ref), atol=1e-12)


def test_principal_angles_ascending_with_multiple_small_angles() -> None:
    """多个小角度同时存在时，混合算法的输出仍须升序（svd 返回降序，需排序）。"""
    rng = np.random.default_rng(21)
    A1 = np.linalg.qr(rng.normal(size=(30, 5)))[0]
    A2 = np.linalg.qr(A1 + 1e-7 * rng.normal(size=(30, 5)))[0]
    theta = principal_angles(A1, A2)
    assert np.all(np.diff(theta) >= -1e-15)
    assert np.allclose(theta, np.sort(subspace_angles(A1, A2)), atol=1e-12)


def test_principal_angles_count_is_min_dim() -> None:
    A1 = random_orthonormal(50, 3, seed=3)
    A2 = random_orthonormal(50, 9, seed=4)
    assert principal_angles(A1, A2).size == 3


def test_identical_subspace_has_zero_angles() -> None:
    A = random_orthonormal(40, 6, seed=5)
    theta = principal_angles(A, A)
    assert np.all(theta < 1e-7)
    assert grassmann_distance(A, A) < 1e-7


def test_distance_symmetry_and_triangle_inequality() -> None:
    A = random_orthonormal(80, 6, seed=6)
    B = random_orthonormal(80, 6, seed=7)
    C = random_orthonormal(80, 6, seed=8)
    assert grassmann_distance(A, B) == pytest.approx(grassmann_distance(B, A))
    assert grassmann_distance(A, C) <= grassmann_distance(A, B) + grassmann_distance(B, C) + 1e-12


def test_distance_invariant_under_ambient_rotation() -> None:
    """子空间距离是内蕴量：环境空间的正交变换不应改变它。"""
    A = random_orthonormal(60, 5, seed=9)
    B = random_orthonormal(60, 5, seed=10)
    Q = np.linalg.qr(np.random.default_rng(11).normal(size=(60, 60)))[0]
    assert grassmann_distance(Q @ A, Q @ B) == pytest.approx(grassmann_distance(A, B), abs=1e-10)


def test_distance_bounded_by_sqrt_k_pi_2() -> None:
    k = 4
    A = random_orthonormal(50, k, seed=12)
    B = random_orthonormal(50, k, seed=13)
    assert grassmann_distance(A, B) <= np.sqrt(k) * np.pi / 2 + 1e-12


def test_distance_matches_scipy_subspace_angles() -> None:
    """维数不同时，距离应等于 scipy 主角度的 ℓ₂ 范数。"""
    A1 = random_orthonormal(50, 4, seed=30)
    A2 = random_orthonormal(50, 6, seed=31)
    ref = float(np.linalg.norm(subspace_angles(A1, A2)))
    assert grassmann_distance(A1, A2) == pytest.approx(ref, abs=1e-12)


def test_asin_branch_improves_small_angle_accuracy() -> None:
    """小角度下 arcsin 分支应显著优于纯 arccos（论文引入混合策略的动机）。"""
    n, k = 40, 4
    theta_true = 1e-7
    A1 = np.eye(n)[:, :k].copy()
    A2 = A1.copy()
    A2[:, 0] = np.cos(theta_true) * np.eye(n)[0] + np.sin(theta_true) * np.eye(n)[k]

    # 主角度升序：被扰动的方向对应最大的那个角，即末位
    hybrid = principal_angles(A1, A2)[-1]
    acos_only = principal_angles(A1, A2, method="acos")[-1]
    assert abs(hybrid - theta_true) < abs(acos_only - theta_true)
    assert hybrid == pytest.approx(theta_true, rel=1e-6)


def test_truncation_rank_monotone_in_p() -> None:
    rng = np.random.default_rng(14)
    M = rng.random((64, 64)) * (rng.random((64, 64)) > 0.4)
    ranks = [truncation_rank(np.linalg.svd(M, compute_uv=False), p) for p in (0.1, 0.3, 0.5, 0.7, 0.9)]
    assert ranks == sorted(ranks)


def test_two_truncation_rules_differ_by_at_most_one() -> None:
    """ge 规则（≥ p·总和）与 le 规则（≤ p·总和）通常只差 1 维。"""
    rng = np.random.default_rng(15)
    for _ in range(20):
        M = rng.random((32, 32))
        s = np.linalg.svd(M, compute_uv=False)
        for p in (0.2, 0.5, 0.71, 0.9):
            a = truncation_rank(s, p, "ge")
            b = truncation_rank(s, p, "le")
            assert abs(a - b) <= 1


def test_basis_is_orthonormal_and_retains_energy() -> None:
    rng = np.random.default_rng(16)
    M = rng.random((48, 48))
    A, info = subspace_basis(M, 0.6, return_info=True)
    assert np.allclose(A.T @ A, np.eye(A.shape[1]), atol=1e-12)
    assert info.retained_ratio >= 0.6
    assert info.rank == A.shape[1]


def test_scale_invariance() -> None:
    """列空间对正数倍缩放不变，故距离也不变（计数矩阵 vs 频率矩阵等价）。"""
    rng = np.random.default_rng(17)
    M1 = rng.poisson(3.0, size=(32, 32)).astype(float)
    M2 = rng.poisson(3.0, size=(32, 32)).astype(float)
    A1 = subspace_basis(M1, 0.7)
    A2 = subspace_basis(M2, 0.7)
    B1 = subspace_basis(M1 * 7.0, 0.7)
    B2 = subspace_basis(M2 * 7.0, 0.7)
    assert grassmann_distance(A1, A2) == pytest.approx(grassmann_distance(B1, B2), abs=1e-10)


def test_invalid_parameters_raise() -> None:
    M = np.ones((8, 8))
    with pytest.raises(ValueError):
        subspace_basis(M, 0.0)
    with pytest.raises(ValueError):
        subspace_basis(M, 1.0)
    with pytest.raises(ValueError):
        subspace_basis(np.zeros((8, 8)), 0.5)
    with pytest.raises(ValueError):
        subspace_basis(M, 0.5, rule="unknown")


def test_default_eps_is_cos_1e_minus_5() -> None:
    assert DEFAULT_EPS == pytest.approx(np.cos(1e-5))
