"""成对距离矩阵的计算（含并行与进度输出）。"""

from __future__ import annotations

import sys
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np

from .fcgr import fcgr_matrix
from .grassmann import (
    DEFAULT_EPS,
    DEFAULT_P,
    basis_from_svd,
    grassmann_distance,
    subspace_basis,
    svd_factors,
)

__all__ = [
    "build_svd_cache",
    "build_bases",
    "distance_matrix_from_bases",
    "distance_matrix",
]


def build_svd_cache(
    sequences: Sequence[str],
    n: int,
    *,
    alphabet: int = 3,
    mode: str = "exact",
    normalize: bool = False,
    verbose: bool = False,
) -> List[Tuple[np.ndarray, np.ndarray]]:
    """对每条序列计算 FCGR 的 SVD 因子 ``(U, s)``，供不同 p 值复用。

    扫描截断参数 p 时无需重复做 FCGR 与 SVD，可显著提速。
    """
    cache: List[Tuple[np.ndarray, np.ndarray]] = []
    t0 = time.perf_counter()
    for i, seq in enumerate(sequences):
        mat = fcgr_matrix(seq, n, alphabet, normalize=normalize, mode=mode)
        cache.append(svd_factors(mat))
        if verbose and (i + 1) % 10 == 0:
            print(f"  [FCGR+SVD] {i + 1}/{len(sequences)}", file=sys.stderr)
    if verbose:
        dt = time.perf_counter() - t0
        print(f"  [FCGR+SVD] 完成 {len(sequences)} 条，用时 {dt:.3f}s", file=sys.stderr)
    return cache


def build_bases(
    sequences: Sequence[str],
    n: int,
    p: float = DEFAULT_P,
    *,
    alphabet: int = 3,
    rule: str = "ge",
    mode: str = "exact",
    normalize: bool = False,
    verbose: bool = False,
    cache: Optional[Sequence[Tuple[np.ndarray, np.ndarray]]] = None,
) -> List[np.ndarray]:
    """对每条序列计算 FCGR 并截断，返回其列空间的正交基列表。

    Parameters
    ----------
    cache : sequence of (U, s), optional
        预先算好的 SVD 因子（见 :func:`build_svd_cache`）；给出时跳过 FCGR 与 SVD。

    Returns
    -------
    list of ndarray
        ``bases[i]`` 形状为 ``(2**n, k_i)``，列正交归一。
    """
    if cache is not None:
        if len(cache) != len(sequences):
            raise ValueError("cache 长度与序列数不一致")
        return [basis_from_svd(U, s, p, rule) for U, s in cache]

    bases: List[np.ndarray] = []
    t0 = time.perf_counter()
    for i, seq in enumerate(sequences):
        mat = fcgr_matrix(seq, n, alphabet, normalize=normalize, mode=mode)
        bases.append(subspace_basis(mat, p, rule=rule))
        if verbose and (i + 1) % 10 == 0:
            print(f"  [FCGR+SVD] {i + 1}/{len(sequences)}", file=sys.stderr)
    if verbose:
        dt = time.perf_counter() - t0
        ranks = [b.shape[1] for b in bases]
        print(
            f"  [FCGR+SVD] 完成 {len(sequences)} 条，用时 {dt:.3f}s，"
            f"截断维数 min/median/max = {min(ranks)}/{int(np.median(ranks))}/{max(ranks)}",
            file=sys.stderr,
        )
    return bases


def distance_matrix_from_bases(
    bases: Sequence[np.ndarray],
    *,
    eps: float = DEFAULT_EPS,
    method: str = "hybrid",
    n_jobs: int = 1,
    verbose: bool = False,
    progress: Optional[Callable[[int, int], None]] = None,
) -> np.ndarray:
    """由子空间基列表计算对称距离矩阵（对角线为 0）。"""
    m = len(bases)
    D = np.zeros((m, m), dtype=np.float64)
    if m < 2:
        return D

    def _row(i: int) -> None:
        Ai = bases[i]
        for j in range(i + 1, m):
            d = grassmann_distance(Ai, bases[j], eps=eps, method=method)
            D[i, j] = d
            D[j, i] = d

    t0 = time.perf_counter()
    if n_jobs is None or n_jobs <= 1:
        for i in range(m - 1):
            _row(i)
            if verbose and (i + 1) % max(1, m // 10) == 0:
                print(f"  [距离] 行 {i + 1}/{m - 1}", file=sys.stderr)
            if progress is not None:
                progress(i + 1, m - 1)
    else:
        # LAPACK 的 SVD 会释放 GIL，因此线程级并行即可获得加速
        with ThreadPoolExecutor(max_workers=n_jobs) as pool:
            list(pool.map(_row, range(m - 1)))
        if progress is not None:
            progress(m - 1, m - 1)
    if verbose:
        print(f"  [距离] 完成 {m * (m - 1) // 2} 对，用时 {time.perf_counter() - t0:.3f}s", file=sys.stderr)
    return D


def distance_matrix(
    sequences: Sequence[str],
    n: int,
    p: float = DEFAULT_P,
    *,
    alphabet: int = 3,
    rule: str = "ge",
    mode: str = "exact",
    normalize: bool = False,
    eps: float = DEFAULT_EPS,
    method: str = "hybrid",
    n_jobs: int = 1,
    verbose: bool = False,
) -> np.ndarray:
    """一步计算序列集合的广义 Grassmann 距离矩阵。

    Parameters
    ----------
    sequences : sequence of str
        DNA 序列集合。
    n : int
        FCGR 分辨率参数（``2**n`` 网格）。
    p : float
        奇异值累积贡献率截断阈值。
    n_jobs : int
        并行线程数，``<=1`` 表示串行。

    Returns
    -------
    ndarray
        形状 ``(m, m)`` 的对称距离矩阵。
    """
    bases = build_bases(
        sequences, n, p, alphabet=alphabet, rule=rule, mode=mode, normalize=normalize, verbose=verbose
    )
    return distance_matrix_from_bases(
        bases, eps=eps, method=method, n_jobs=n_jobs, verbose=verbose
    )
