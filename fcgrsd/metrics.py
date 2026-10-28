"""树拓扑距离（Robinson-Foulds）与经典多维标度（MDS）。"""

from __future__ import annotations

from typing import Tuple

import numpy as np

from .tree import Tree

__all__ = ["robinson_foulds", "classical_mds"]


def robinson_foulds(
    tree_a: Tree,
    tree_b: Tree,
    *,
    normalize: bool = False,
    rooted: bool = False,
) -> float:
    """Robinson-Foulds 距离（拓扑差异）。

    Parameters
    ----------
    tree_a, tree_b : Tree
        待比较的两棵树，叶子名集合必须一致。
    normalize : bool
        为 ``True`` 时返回 ``RF / (|S_a| + |S_b|)``，取值 ``[0, 1]``。
    rooted : bool
        为 ``True`` 时按有根树的 clade 集合比较，否则按无根树的二分划比较。

    Returns
    -------
    float
        RF 距离（整数计数；``normalize=True`` 时为归一化值）。
    """
    if rooted:
        sa = tree_a.clusters()
        sb = tree_b.clusters()
    else:
        sa = tree_a.splits()
        sb = tree_b.splits()

    la = {n.name for n in tree_a.root.leaves()}
    lb = {n.name for n in tree_b.root.leaves()}
    if len(la) != tree_a.n_leaves or len(lb) != tree_b.n_leaves:
        raise ValueError("树中存在重复叶名，无法计算 RF 距离")
    if la != lb:
        only_a = sorted(x for x in la - lb if x)
        only_b = sorted(x for x in lb - la if x)
        raise ValueError(
            "两棵树的叶子集合不一致，无法计算 RF 距离："
            f"仅在 A 中 {only_a[:5]}，仅在 B 中 {only_b[:5]}"
        )

    diff = len(sa ^ sb)
    if normalize:
        denom = len(sa) + len(sb)
        return float(diff) / denom if denom else 0.0
    return float(diff)


def classical_mds(
    distances: np.ndarray,
    n_components: int = 3,
    *,
    positive: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """经典多维标度（Torgerson / PCoA）。

    论文用 MDS 把 Grassmann 距离矩阵嵌入三维欧氏空间做点云可视化
    （Figure 4–6），此处给出同一算法的实现。

    Parameters
    ----------
    distances : ndarray
        形状 ``(m, m)`` 的对称距离矩阵。
    n_components : int
        目标维数。
    positive : bool
        为 ``True`` 时把负特征值截断为 0（保证取平方根合法）。

    Returns
    -------
    coords : ndarray
        形状 ``(m, min(n_components, m))`` 的坐标，按主坐标方差降序排列。
    eigenvalues : ndarray
        形状 ``(min(n_components, m),)`` 的特征值。
    """
    D = np.asarray(distances, dtype=np.float64)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError("distances 必须是方阵")
    if not np.all(np.isfinite(D)):
        raise ValueError("distances 包含 NaN 或 Inf")
    if not np.allclose(D, D.T, atol=1e-12):
        D = 0.5 * (D + D.T)
    m = D.shape[0]
    if int(n_components) < 1:
        raise ValueError("n_components 必须为正整数")
    k = min(int(n_components), m)
    D2 = D**2
    J = np.eye(m) - np.ones((m, m)) / m
    B = -0.5 * J @ D2 @ J
    B = 0.5 * (B + B.T)
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1][:k]
    vals = vals[order]
    vecs = vecs[:, order]
    if positive:
        vals = np.clip(vals, 0.0, None)
    coords = vecs * np.sqrt(vals)
    return coords, vals
