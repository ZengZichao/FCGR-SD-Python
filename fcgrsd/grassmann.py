"""Grassmann-manifold building blocks: truncated subspaces, principal angles, distance.

Pipeline (paper supplementary §1.2.2–1.2.3 and Algorithm 3)
----------------------------------------------------------
**Step 1: truncated SVD.** For an FCGR matrix :math:`M` (size
:math:`2^n\\times2^n`),

.. math::  M = P D Q^{\\mathsf T},

keep the first k left singular vectors

.. math::  A = P_{:,1:k},

where k is selected by the cumulative-energy rule implemented in
:func:`truncation_rank`. ``rule="ge"`` (default) takes the smallest k whose
cumulative sum reaches :math:`p\\sum_j d_j`; ``rule="le"`` takes the largest k
whose cumulative sum does not exceed it. The two differ by at most one
dimension; ``"ge"`` reproduces the published numbers.

**Step 2: principal angles.** For orthonormal bases
:math:`A_1\\in\\mathbb R^{N\\times k_1}` and :math:`A_2\\in\\mathbb R^{N\\times k_2}`,

.. math::  A_1^{\\mathsf T}A_2 = U\\Sigma V^{\\mathsf T},
           \\qquad \\theta_i = \\arccos\\sigma_i .

**Step 3: generalized Grassmann distance.**

.. math::  d(A_1,A_2) = \\Big(\\sum_{i=1}^{\\min(k_1,k_2)} \\theta_i^2\\Big)^{1/2}.

Numerical detail: :math:`\\arccos` loses precision as :math:`\\sigma_i\\to1`
(relative error about :math:`\\epsilon/\\sqrt{1-\\sigma_i^2}`). When
:math:`\\sigma_i \\ge \\cos(10^{-5})` (i.e. :math:`\\theta_i \\lesssim 10^{-5}`
radians) the sine form is used instead:

.. math::  \\mu_i = \\sigma_i\\big((I-A_1A_1^{\\mathsf T})A_2 v_i\\big) = \\sin\\theta_i,
           \\qquad \\theta_i = \\arcsin\\mu_i .
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Tuple

import numpy as np

__all__ = [
    "DEFAULT_EPS",
    "TruncationInfo",
    "subspace_basis",
    "truncated_basis",
    "truncation_rank",
    "basis_from_svd",
    "svd_factors",
    "principal_angles",
    "grassmann_distance",
    "DEFAULT_P",
]

#: Cosine threshold for the small-angle branch: ``cos(1e-5)``.
DEFAULT_EPS: float = float(np.cos(1e-5))

#: Optimal truncation threshold reported by the paper for 44 SARS-CoV-2 genomes.
DEFAULT_P: float = 0.75

TruncationRule = Literal["ge", "le"]


@dataclass(frozen=True)
class TruncationInfo:
    """截断 SVD 的元信息。"""

    rank: int
    """截断后保留的维数 k。"""

    total_singular_sum: float
    """全部奇异值之和。"""

    retained_ratio: float
    """前 k 个奇异值之和占总和的比例。"""

    singular_values: np.ndarray
    """降序排列的全部奇异值。"""


def _as_matrix(M: np.ndarray) -> np.ndarray:
    arr = np.asarray(M)
    if arr.ndim != 2:
        raise ValueError(f"期望二维矩阵，收到 ndim={arr.ndim}")
    return arr


def truncation_rank(s: np.ndarray, p: float, rule: TruncationRule = "ge") -> int:
    """Determine the truncation rank k from the cumulative singular-value energy.

    ``rule="ge"``: smallest k with :math:`\\sum_{j\\le k} s_j \\ge p\\sum_j s_j`
    (default; reproduces the published results).

    ``rule="le"``: largest k with :math:`\\sum_{j\\le k} s_j \\le p\\sum_j s_j`
    (the literal wording of the paper). If no k satisfies it, the rank falls
    back to 1.
    """
    if not 0.0 < p < 1.0:
        raise ValueError(f"p must lie in (0, 1), got {p!r}")
    s = np.asarray(s, dtype=np.float64)
    total = float(s.sum())
    if total <= 0.0:
        raise ValueError("all singular values are zero; cannot define a column space")
    target = p * total
    cum = np.cumsum(s)
    if rule == "ge":
        k = int(np.searchsorted(cum, target, side="left")) + 1
    elif rule == "le":
        k = int(np.searchsorted(cum, target, side="right"))
    else:
        raise ValueError("rule must be 'ge' or 'le'")
    return max(1, min(k, s.size))


def basis_from_svd(
    U: np.ndarray,
    s: np.ndarray,
    p: float = DEFAULT_P,
    rule: TruncationRule = "ge",
) -> np.ndarray:
    """由预先算好的 SVD 因子直接截取正交基（供 p 值扫描复用，避免重复 SVD）。"""
    k = truncation_rank(s, p, rule)
    return np.ascontiguousarray(U[:, :k])


def svd_factors(M: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """返回 FCGR 的左奇异向量与奇异值（降序）。"""
    arr = _as_matrix(M)
    U, s, _ = np.linalg.svd(arr, full_matrices=False)
    return U, s


def subspace_basis(
    M: np.ndarray,
    p: float = DEFAULT_P,
    *,
    rule: TruncationRule = "ge",
    return_info: bool = False,
):
    """按奇异值累积贡献率截断，返回列空间的标准正交基。

    Parameters
    ----------
    M : ndarray
        FCGR 矩阵（计数或频率均可；正数倍缩放不改变列空间）。
    p : float
        累积贡献率阈值，取值 ``(0, 1)``。论文推荐 0.71（30 条 HA 基因，
        n=6）、0.75（44 条 SARS-CoV-2 基因组，n=8）、0.43（3283 条冠状病毒）。
    rule : {"ge", "le"}
        Truncation rule; the two differ by at most one dimension:

        - ``"ge"`` (default): smallest k with
          :math:`\\sum_{j\\le k}d_j \\ge p\\sum_j d_j`;
        - ``"le"``: largest k with
          :math:`\\sum_{j\\le k}d_j \\le p\\sum_j d_j` (the literal wording of
          the paper's supplementary §1.2.3).
    return_info : bool
        是否同时返回 :class:`TruncationInfo`。

    Returns
    -------
    A : ndarray
        形状 ``(N, k)``，列正交归一。
    info : TruncationInfo, optional
        仅当 ``return_info=True`` 时返回。
    """
    arr = _as_matrix(M)
    U, s, _ = np.linalg.svd(arr, full_matrices=False)
    total = float(s.sum())
    k = truncation_rank(s, p, rule)
    A = np.ascontiguousarray(U[:, :k])
    if not return_info:
        return A
    cum = np.cumsum(s)
    info = TruncationInfo(
        rank=k,
        total_singular_sum=total,
        retained_ratio=float(cum[k - 1] / total),
        singular_values=s,
    )
    return A, info


def truncated_basis(M: np.ndarray, p: float = DEFAULT_P, **kwargs):
    """``subspace_basis`` 的别名，语义上更贴近论文中的"截断特征矩阵"。"""
    return subspace_basis(M, p, **kwargs)


def principal_angles(
    A1: np.ndarray,
    A2: np.ndarray,
    *,
    eps: float = DEFAULT_EPS,
    method: Literal["hybrid", "acos"] = "hybrid",
) -> np.ndarray:
    """计算两个子空间之间的主角度（弧度，升序）。

    Parameters
    ----------
    A1, A2 : ndarray
        形状分别为 ``(N, k1)``、``(N, k2)`` 的列正交归一矩阵。
    eps : float
        余弦阈值，默认 :math:`\\cos(10^{-5})`。当 :math:`\\sigma_i \\ge \\epsilon`
        （即角度足够小）时改用 :math:`\\arcsin` 分支以保证精度。
    method : {"hybrid", "acos"}
        ``hybrid``：论文 Algorithm 3 的混合策略；
        ``acos``：全程 :math:`\\arccos`（更慢但可用于对照）。

    Returns
    -------
    ndarray
        长度为 ``min(k1, k2)`` 的主角度数组，已裁剪到 ``[0, pi/2]``。
    """
    A1 = _as_matrix(A1)
    A2 = _as_matrix(A2)
    if A1.shape[0] != A2.shape[0]:
        raise ValueError(
            f"两个基的行数（环境空间维数）必须一致，收到 {A1.shape} 与 {A2.shape}"
        )

    inner = A1.T @ A2
    _, s, Vt = np.linalg.svd(inner, full_matrices=False)
    s = np.clip(s, 0.0, 1.0)
    p = s.size
    if p == 0:
        return np.zeros(0)

    if method == "acos":
        return np.arccos(s)
    if method != "hybrid":
        raise ValueError("method 只能是 'hybrid' 或 'acos'")

    # 若最大奇异值都低于阈值，说明最小主角度已超过 1e-5，arccos 足够精确
    if s[0] < eps:
        return np.arccos(s)

    below = np.flatnonzero(s < eps)
    n_small = int(below[0]) if below.size else p  # 前 n_small 个角度属于"小角度"

    theta = np.empty(p, dtype=np.float64)
    theta[n_small:] = np.arccos(s[n_small:])
    if n_small:
        # 正弦形式：mu_i = sigma_i((I - A1 A1^T) A2 v_i)
        Vcos = A2 @ Vt.T
        R = Vcos[:, :n_small]
        B = R - A1 @ (A1.T @ R)
        # svd 返回降序，而主角度约定为升序，故排序后再取 arcsin
        mu = np.sort(np.linalg.svd(B, compute_uv=False)[:n_small])
        theta[:n_small] = np.arcsin(np.clip(mu, 0.0, 1.0))
    return np.clip(theta, 0.0, np.pi / 2)


def grassmann_distance(
    A1: np.ndarray,
    A2: np.ndarray,
    *,
    eps: float = DEFAULT_EPS,
    method: Literal["hybrid", "acos"] = "hybrid",
) -> float:
    """两个子空间之间的广义 Grassmann 距离 :math:`\\|\\theta\\|_2`。

    当两个子空间维数不同时，只取 :math:`\\min(k_1,k_2)` 个主角度——
    这与主角度本身的定义一致（个数等于两个维数的较小者）。
    """
    theta = principal_angles(A1, A2, eps=eps, method=method)
    return float(np.sqrt(np.dot(theta, theta)))
