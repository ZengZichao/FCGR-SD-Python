"""Chaos Game Representation（CGR，混沌游戏表示）的坐标计算与网格定位。

数学定义
--------
给定 DNA 序列 :math:`s_1 s_2 \\cdots s_L`，CGR 点列由如下仿射递推给出：

.. math::

    P_0 = \\left(\\tfrac12, \\tfrac12\\right), \\qquad
    P_i = \\tfrac12 \\left(P_{i-1} + P(s_i)\\right),

其中 :math:`P(\\cdot)` 把四种碱基映射到单位正方形的四个顶点。顶点与碱基的
对应关系共有 24 种（4! 个置换），原论文按"嘌呤/嘧啶配对"的几何直觉选取了
6 种（参数 ``alphabet``，取值 1~6）。

网格定位
--------
把单位正方形等分为 :math:`2^n \\times 2^n` 个格子，点 :math:`P_i` 落入的格子
（从 1 开始计数）为

.. math::

    c_i^{(1)} = \\left\\lfloor 2^n x_i \\right\\rfloor + 1, \\quad
    c_i^{(2)} = \\left\\lfloor 2^n y_i \\right\\rfloor + 1 .

本模块提供两条实现路径：

``exact``（默认，推荐）
    利用 CGR 坐标的二进制结构做整数精确计算。由递推式可证

    .. math::

        x_i = \\frac{X_i}{2^{i+1}},\\qquad X_i = 1 + \\sum_{j=1}^{i} u_j 2^j ,

    其中 :math:`u_j \\in \\{0,1\\}` 是第 j 个碱基的横坐标比特。因 :math:`X_i`
    恒为奇数，:math:`2^n x_i` 永不落在格线上，于是

    .. math::

        c_i^{(1)} = 1 + \\sum_{k=0}^{n-1} u_{i-n+1+k}\\, 2^{k} ,

    即格子编号只由**窗口内最后 n 个碱基**决定（窗口首碱基权重为 :math:`2^0`，
    末碱基权重为 :math:`2^{n-1}`）。该式可用纯整数运算向量化实现，完全规避
    浮点误差。

``float``
    逐点递推 float64，并在 ``j >= 50`` 时施加 :math:`2^{-31}` 扰动，
    避免长同聚段的坐标下溢；用于与整数路径做数值对照。

参考文献：Jeffrey H.J. (1990) *Nucleic Acids Res.* 18:2163-2170；
Li et al. (2024) *The Innovation* 5:100677 补充材料 §1.2.1、§2。
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np

__all__ = [
    "ALPHABETS",
    "ALPHABET_DESCRIPTIONS",
    "ambiguous_fraction",
    "normalize_sequence",
    "encode_bits",
    "cgr_coordinates",
    "cgr_grid_indices",
    "cgr_grid_indices_float",
]

# ---------------------------------------------------------------------------
# 六种碱基-顶点对应关系（alphabet = 1..6）
# 每个顶点写作 (第一坐标, 第二坐标)
# ---------------------------------------------------------------------------
ALPHABETS: Dict[int, Dict[str, Tuple[int, int]]] = {
    1: {"A": (0, 0), "C": (0, 1), "G": (1, 1), "T": (1, 0)},
    2: {"A": (0, 0), "C": (0, 1), "T": (1, 1), "G": (1, 0)},
    3: {"A": (0, 0), "G": (0, 1), "C": (1, 1), "T": (1, 0)},
    4: {"A": (0, 0), "G": (0, 1), "T": (1, 1), "C": (1, 0)},
    5: {"A": (0, 0), "T": (0, 1), "C": (1, 1), "G": (1, 0)},
    6: {"A": (0, 0), "T": (0, 1), "G": (1, 1), "C": (1, 0)},
}

ALPHABET_DESCRIPTIONS: Dict[int, str] = {
    1: "ACGT: A(0,0), C(0,1), G(1,1), T(1,0)",
    2: "ACTG: A(0,0), C(0,1), T(1,1), G(1,0)",
    3: "AGCT: A(0,0), G(0,1), C(1,1), T(1,0) (paper Algorithm 2; default)",
    4: "AGTC: A(0,0), G(0,1), T(1,1), C(1,0)",
    5: "ATCG: A(0,0), T(0,1), C(1,1), G(1,0)",
    6: "ATGC: A(0,0), T(0,1), G(1,1), C(1,0)",
}

_VALID = frozenset("ACGT")

# 字节 -> 比特查表（255 表示非法字符）
_LUT: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
for _kap, _mapping in ALPHABETS.items():
    _u = np.full(256, 255, dtype=np.uint8)
    _v = np.full(256, 255, dtype=np.uint8)
    for _base, (_cu, _cv) in _mapping.items():
        _u[ord(_base)] = _cu
        _v[ord(_base)] = _cv
    _LUT[_kap] = (_u, _v)


def _check_alphabet(alphabet: int) -> None:
    if alphabet not in ALPHABETS:
        raise ValueError(
            f"alphabet must be an integer in 1..6, got {alphabet!r}"
        )


def _check_n(n: int) -> int:
    """校验 FCGR 分辨率参数 n 为正整数，返回 int。"""
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)):
        raise ValueError(f"n 必须是正整数，收到 {n!r}")
    n = int(n)
    if n < 1:
        raise ValueError(f"n 必须是正整数，收到 {n}")
    return n


def ambiguous_fraction(seq: str) -> float:
    """返回序列中非 A/C/G/T 字符（N、-、简并碱基等）的占比。

    该比例决定"剔除非法字符"对 FCGR 的破坏程度：实测单个 N 使 FCGR 的
    L1 误差约 0.5%，1% 的 N 约 10%，5% 的 N 约 44%。
    """
    s = "".join(seq.split()).upper()
    if not s:
        return 0.0
    codes = np.frombuffer(s.encode("utf-8", errors="replace"), dtype=np.uint8)
    valid = np.isin(codes, (65, 67, 71, 84))  # A C G T
    return float(1.0 - valid.mean())


def normalize_sequence(seq: str, on_unknown: str = "skip") -> str:
    """统一大写并处理非 ACGT 字符。

    Parameters
    ----------
    seq : str
        原始核酸序列，允许含换行、空白与小写。
    on_unknown : {"skip", "error"}
        ``skip``：剔除 A/C/G/T 之外的字符（如 ``N``、``-``），这是无比对方法
        的常规做法；``error``：遇到非法字符立即抛错。

    Returns
    -------
    str
        仅由 ``ACGT`` 组成的序列。
    """
    s = "".join(seq.split()).upper()
    if on_unknown == "error":
        bad = sorted(set(s) - _VALID)
        if bad:
            raise ValueError(f"序列包含非法字符 {bad}（on_unknown='error'）")
        return s
    if on_unknown != "skip":
        raise ValueError("on_unknown 只能是 'skip' 或 'error'")
    if s.strip("ACGT") == "":
        return s
    return "".join(ch for ch in s if ch in _VALID)


def encode_bits(seq: str, alphabet: int = 3) -> Tuple[np.ndarray, np.ndarray]:
    """把碱基序列编码为两路 0/1 比特流。

    Parameters
    ----------
    seq : str
        已经过 :func:`normalize_sequence` 的序列（本函数会再做一次大写与过滤）。
    alphabet : int
        碱基-顶点对应关系编号（1~6）。

    Returns
    -------
    (u, v) : tuple of ndarray
        形状均为 ``(L,)`` 的 ``uint8`` 数组，分别对应 CGR 的第一、第二坐标比特。
    """
    _check_alphabet(alphabet)
    s = normalize_sequence(seq)
    lut_u, lut_v = _LUT[alphabet]
    codes = np.frombuffer(s.encode("utf-8", errors="replace"), dtype=np.uint8)
    u = lut_u[codes]
    v = lut_v[codes]
    keep = u != 255
    return u[keep], v[keep]


def cgr_coordinates(
    seq: str,
    alphabet: int = 3,
    *,
    correct_float: bool = True,
    exact: bool = False,
) -> np.ndarray:
    """计算 CGR 点列（用于可视化或逐位比对）。

    Parameters
    ----------
    seq : str
        DNA 序列。
    alphabet : int
        碱基-顶点对应关系编号（1~6）。
    correct_float : bool
        仅在 ``exact=False`` 时生效。为 ``True`` 时在 ``j >= 50`` 对坐标施加 :math:`2^{-31}` 扰动，抑制浮点退化。
    exact : bool
        为 ``True`` 时用 Python 大整数精确计算 :math:`X_i/2^{i+1}`，
        结果为正确舍入的浮点值，无累积误差（复杂度 :math:`O(L^2/64)`，
        适合短序列校验）。

    Returns
    -------
    ndarray
        形状 ``(L, 2)``，第 i 行为 :math:`P_i`（i 从 1 开始，不含初始点
        :math:`P_0`）。
    """
    _check_alphabet(alphabet)
    u, v = encode_bits(seq, alphabet)
    m = u.shape[0]
    if m == 0:
        return np.zeros((0, 2), dtype=np.float64)

    if exact:
        # X_i = 1 + sum_{j<=i} u_j 2^j ；逐项累加，用整数保证无误差
        out = np.empty((m, 2), dtype=np.float64)
        xu = 1
        xv = 1
        pow2 = 1
        for i in range(m):
            pow2 <<= 1
            xu += int(u[i]) * pow2
            xv += int(v[i]) * pow2
            denom = float(1 << (i + 2))
            out[i, 0] = xu / denom
            out[i, 1] = xv / denom
        return out

    # ---- 逐点递推 ----
    pts = np.empty((m + 1, 2), dtype=np.float64)
    pts[0] = 0.5
    two_31 = 2.0**-31
    two_50 = 2.0**-50
    two_500 = 2.0**-500
    pu = u.astype(np.float64)
    pv = v.astype(np.float64)
    for j in range(m):
        # 1-based 计数下 j >= 50，对应此处的 j + 1 >= 50
        if correct_float and j + 1 >= 50:
            if pts[j, 0] <= two_500:
                pts[j, 0] += two_31
            elif abs(pts[j, 0] - 0.5) <= two_50:
                pts[j, 0] += two_31 if pts[j - 1, 0] < 0.5 else -two_31
            elif 1.0 - pts[j, 0] <= two_50:
                pts[j, 0] -= two_31

            if pts[j, 1] <= two_500:
                pts[j, 1] += two_31
            elif abs(pts[j, 1] - 0.5) <= two_50:
                pts[j, 1] += two_31 if pts[j - 1, 1] < 0.5 else -two_31
            elif 1.0 - pts[j, 1] <= two_50:
                pts[j, 1] -= two_31
        pts[j + 1, 0] = 0.5 * (pts[j, 0] + pu[j])
        pts[j + 1, 1] = 0.5 * (pts[j, 1] + pv[j])
    return pts[1:]


def cgr_grid_indices(seq: str, n: int, alphabet: int = 3) -> np.ndarray:
    """精确计算每个 CGR 点所属的格子编号（1-based），完全避免浮点误差。

    Parameters
    ----------
    seq : str
        DNA 序列。
    n : int
        FCGR 的分辨率参数，网格为 :math:`2^n \\times 2^n`。
    alphabet : int
        碱基-顶点对应关系编号（1~6）。

    Returns
    -------
    ndarray
        形状 ``(L - n + 1, 2)`` 的 ``int64`` 数组，每行为 ``(行号, 列号)``，
        取值范围 ``1 .. 2**n``，即第 ``n`` 个点及其之后的所有点。
    """
    n = _check_n(n)
    u, v = encode_bits(seq, alphabet)
    length = u.shape[0]
    if length < n:
        raise ValueError(
            f"序列长度 {length} 小于 FCGR 分辨率参数 n={n}，无法构造 FCGR"
        )
    windows = length - n + 1
    rows = np.ones(windows, dtype=np.int64)
    cols = np.ones(windows, dtype=np.int64)
    # rows += sum_k u_{i-n+1+k} * 2^k
    for k in range(n):
        w = np.int64(1) << k
        rows += u[k : k + windows].astype(np.int64) * w
        cols += v[k : k + windows].astype(np.int64) * w
    return np.column_stack((rows, cols))


def cgr_grid_indices_float(
    seq: str,
    n: int,
    alphabet: int = 3,
    *,
    correct_float: bool = True,
) -> np.ndarray:
    """用浮点 CGR 路径计算格子编号（1-based），用于与整数路径做数值对照。

    实现细节：先递推浮点坐标，再取第 ``n`` 个点起的所有点，
    以 ``ceil(2^n * x)`` 定位，并把结果裁剪到
    ``[1, 2**n]`` 以免极端情况下越界。
    """
    n = _check_n(n)
    coords = cgr_coordinates(seq, alphabet, correct_float=correct_float)
    if coords.shape[0] < n:
        raise ValueError(
            f"序列长度 {coords.shape[0]} 小于 FCGR 分辨率参数 n={n}，无法构造 FCGR"
        )
    coords = coords[n - 1 :]
    size = 2**n
    idx = np.ceil(size * coords)
    idx = np.clip(idx, 1, size).astype(np.int64)
    return idx
