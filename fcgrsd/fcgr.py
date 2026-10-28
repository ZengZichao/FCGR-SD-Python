"""FCGR (Frequency matrix of Chaos Game Representation).

A sequence of length L is mapped through CGR to a :math:`2^n \\times 2^n`
count matrix. ``fcgr[i, j]`` is the number of CGR points falling in grid cell
``(i, j)`` (0-based).

Two grid-location paths are available:

``exact`` (default)
    Integer closed form derived from the binary structure of the CGR
    coordinates (see :mod:`fcgrsd.cgr`). No floating-point rounding.

``float``
    The CGR recurrence in float64 followed by ``ceil(2^n x)``, clipped to the
    grid. Useful as a numerical cross-check of the integer path.
"""

from __future__ import annotations

from typing import List, Sequence

import numpy as np

from .cgr import (
    _check_n,
    ambiguous_fraction,
    cgr_grid_indices,
    cgr_grid_indices_float,
    normalize_sequence,
)

__all__ = ["fcgr_matrix", "fcgr_matrices", "fcgr_total"]

_MAX_RECOMMENDED_N = 11

#: 非 ACGT 字符占比超过该阈值时给出警告
_AMBIGUOUS_WARN_THRESHOLD = 0.01


def _warn_ambiguous(fraction: float) -> None:
    """提示非 ACGT 字符占比过高。

    不做全局去重：Python 的 warnings 机制本身会按 (消息, 位置) 去重，
    而消息里带有具体占比，因此不同严重程度的序列都会各自提示一次。
    早期的全局标志会让第一条之后的序列完全静默，掩盖更严重的输入。
    """
    import warnings

    warnings.warn(
        f"序列中非 ACGT 字符占比 {fraction:.1%}，已被剔除。该比例越高，FCGR 失真越明显"
        "（实测 1% → L1 误差约 10%，5% → 约 44%）。建议先过滤此类序列，"
        "或用 on_unknown='error' 强制报错。",
        UserWarning,
        stacklevel=3,
    )


def fcgr_matrix(
    seq: str,
    n: int,
    alphabet: int = 3,
    *,
    normalize: bool = False,
    mode: str = "exact",
    dtype=None,
) -> np.ndarray:
    """计算单条序列的 FCGR 矩阵。

    Parameters
    ----------
    seq : str
        DNA 序列（大小写不敏感，非 ACGT 字符默认剔除）。
    n : int
        FCGR 分辨率参数，矩阵尺寸为 ``2**n``。
    alphabet : int
        碱基-顶点对应关系编号（1~6），见 :data:`fcgrsd.cgr.ALPHABETS`。
    normalize : bool
        为 ``True`` 时把计数矩阵除以总计数，得到"频率矩阵"（原论文的称谓）。
        注意：对本方法而言归一化不改变列空间，因而不改变任何距离。
    mode : {"exact", "float"}
        ``exact``：整数精确网格定位（默认，推荐）；
        ``float``：浮点 CGR 递推 + ``ceil`` 定位。
    dtype : numpy dtype, optional
        输出数组类型。默认计数为 ``int64``，归一化后为 ``float64``。

    Returns
    -------
    ndarray
        形状 ``(2**n, 2**n)`` 的矩阵。
    """
    n = _check_n(n)
    _warn_large_n(n)
    frac = ambiguous_fraction(seq)
    if frac > _AMBIGUOUS_WARN_THRESHOLD:
        _warn_ambiguous(frac)
    size = 1 << n
    if mode == "exact":
        idx = cgr_grid_indices(seq, n, alphabet)
    elif mode == "float":
        idx = cgr_grid_indices_float(seq, n, alphabet)
    else:
        raise ValueError("mode 只能是 'exact' 或 'float'")

    flat = (idx[:, 0] - 1) * size + (idx[:, 1] - 1)
    counts = np.bincount(flat, minlength=size * size)
    mat = counts.reshape(size, size)

    if normalize:
        total = mat.sum()
        out = mat.astype(np.float64)
        if total > 0:
            out /= total
        return out if dtype is None else out.astype(dtype, copy=False)
    if dtype is not None:
        return mat.astype(dtype, copy=False)
    return mat


def fcgr_matrices(
    sequences: Sequence[str],
    n: int,
    alphabet: int = 3,
    *,
    normalize: bool = False,
    mode: str = "exact",
    dtype=None,
) -> List[np.ndarray]:
    """批量计算 FCGR，返回列表（不堆叠成三维数组以节省内存）。"""
    return [
        fcgr_matrix(s, n, alphabet, normalize=normalize, mode=mode, dtype=dtype)
        for s in sequences
    ]


def fcgr_total(seq: str, n: int) -> int:
    """返回参与计数的 CGR 点数，即 ``len(clean(seq)) - n + 1``。

    ``clean(seq)`` 是剔除 A/C/G/T 之外字符后的序列；因此含 ``N``/``-`` 的
    序列不会高估计数。
    """
    n = _check_n(n)
    return max(0, len(normalize_sequence(seq)) - n + 1)


def _warn_large_n(n: int) -> None:
    import warnings

    if n > _MAX_RECOMMENDED_N:
        warnings.warn(
            f"n={n} 会生成 {1 << n} x {1 << n} 的稠密矩阵，内存与 SVD 开销增长很快；"
            "经验上应取 n ≈ log4(L)，L 为序列长度。",
            RuntimeWarning,
            stacklevel=3,
        )
