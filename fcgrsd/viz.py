"""可视化工具：CGR 点云、FCGR 热图、矩形进化树、MDS 散点。

matplotlib 为可选依赖，仅在调用本模块的函数时导入。
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np

from .cgr import cgr_coordinates
from .tree import Node, Tree

__all__ = ["plot_cgr", "plot_fcgr", "draw_tree", "plot_mds"]


def _plt():
    import matplotlib

    if not matplotlib.get_backend().lower().startswith(("agg", "pdf", "svg", "ps")):
        matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    return plt


def plot_cgr(
    seq: str,
    ax=None,
    *,
    alphabet: int = 3,
    s: float = 0.6,
    cmap: str = "viridis",
    walk: int = 0,
):
    """绘制 CGR 点云。

    Parameters
    ----------
    walk : int
        大于 0 时，额外用折线画出前 ``walk`` 个点的"走笔"轨迹，
        便于观察 CGR 的分形结构（短序列点云本身很稀疏）。
    """
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(4.2, 4.2), dpi=160)
    coords = cgr_coordinates(seq, alphabet)
    if walk and walk > 1:
        k = min(walk, len(coords))
        xs = np.concatenate(([0.5], coords[:k, 0]))
        ys = np.concatenate(([0.5], coords[:k, 1]))
        ax.plot(xs, ys, color="#4a5568", lw=0.5, alpha=0.65, zorder=1)
        ax.scatter(xs[1:], ys[1:], c=np.arange(k), cmap=cmap, s=5, linewidths=0, zorder=2)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], c=np.arange(len(coords)), cmap=cmap, s=s, linewidths=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("CGR" + (f" (first {walk} bases)" if walk else ""), fontsize=10)
    return ax


def plot_fcgr(matrix: np.ndarray, ax=None, *, log: bool = True, cmap: str = "magma"):
    """绘制 FCGR 热图（默认做 log(1+x) 变换以便观察稀疏结构）。

    注意：``fcgr[i, j]`` 的第 i 维是 CGR 的第一坐标（横轴）、第 j 维是第二坐标
    （纵轴）；``imshow`` 默认把第 0 维放在纵轴，因此这里转置后再画，保证与
    :func:`plot_cgr` 的朝向一致。
    """
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(4.2, 4.2), dpi=160)
    mat = np.asarray(matrix).T
    data = np.log1p(mat) if log else mat
    im = ax.imshow(data, origin="lower", cmap=cmap, interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"FCGR {matrix.shape[0]}×{matrix.shape[1]}", fontsize=10)
    return im


def draw_tree(
    tree: Tree,
    ax=None,
    *,
    colors: Optional[Dict[str, str]] = None,
    label_size: float = 6.0,
    line_width: float = 1.0,
    figsize: Optional[tuple] = None,
):
    """绘制矩形（直角）系统发育树。

    Parameters
    ----------
    colors : dict, optional
        叶子名 -> 颜色；用于按类群着色（对应论文 Figure 1/2 的配色）。
    """
    plt = _plt()
    leaves: List[Node] = tree.root.leaves()
    m = len(leaves)
    ypos: Dict[int, float] = {id(leaf): float(i) for i, leaf in enumerate(leaves)}
    xpos: Dict[int, float] = {}

    def assign_x(node: Node, x: float) -> None:
        xpos[id(node)] = x
        for c in node.children:
            assign_x(c, x + (c.length or 0.0))

    def assign_y(node: Node) -> float:
        if node.is_leaf:
            return ypos[id(node)]
        ys = [assign_y(c) for c in node.children]
        ypos[id(node)] = float(np.mean(ys))
        return ypos[id(node)]

    assign_x(tree.root, 0.0)
    assign_y(tree.root)

    if ax is None:
        if figsize is None:
            figsize = (6.2, max(2.6, 0.16 * m + 1.2))
        _, ax = plt.subplots(figsize=figsize, dpi=170)

    def draw(node: Node) -> None:
        x0 = xpos[id(node)]
        if node.is_leaf:
            return
        ys = [ypos[id(c)] for c in node.children]
        ax.plot(
            [x0, x0],
            [min(ys), max(ys)],
            color=_branch_color(node),
            lw=line_width,
            solid_capstyle="butt",
        )
        for c in node.children:
            x1 = xpos[id(c)]
            y1 = ypos[id(c)]
            ax.plot([x0, x1], [y1, y1], color=_branch_color(c), lw=line_width, solid_capstyle="butt")
            draw(c)

    def _branch_color(node: Node) -> str:
        """若某分支下的所有叶子同色，则整枝着色（便于识别单系群）。"""
        if not colors:
            return "#333333"
        cols = {colors.get(leaf.name) for leaf in node.leaves()}
        return cols.pop() if len(cols) == 1 and None not in cols else "#333333"

    draw(tree.root)
    xmax = max(xpos.values())
    xmin = min(xpos.values())
    span = xmax - xmin
    if span <= 0:
        # 所有枝长为 0（例如完全相同的序列）时，给一个单位宽度，避免退化坐标轴
        span = 1.0
        xmax = xmin + span
    for i, leaf in enumerate(leaves):
        col = colors.get(leaf.name, "#222222") if colors else "#222222"
        ax.text(xmax + 0.01 * span, i, leaf.name, va="center", ha="left", fontsize=label_size, color=col)
    ax.set_xlim(xmin - 0.01 * span, xmax + 0.22 * span)
    ax.set_ylim(-1, m)
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("Grassmann distance", fontsize=8)
    ax.tick_params(axis="x", labelsize=7)
    return ax


def plot_mds(
    coords: np.ndarray,
    names: Sequence[str],
    ax=None,
    *,
    colors: Optional[Dict[str, str]] = None,
    dims: int = 3,
    s: float = 26,
    annotate: bool = False,
    label_size: float = 5.5,
):
    """MDS 点云散点图（3D 或 2D）。"""
    plt = _plt()
    if ax is None:
        fig = plt.figure(figsize=(5.6, 4.8), dpi=170)
        ax = fig.add_subplot(111, projection="3d" if dims >= 3 else None)
    cols = [colors.get(n, "#2b6cb0") for n in names] if colors else "#2b6cb0"
    if dims >= 3:
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=cols, s=s, depthshade=False)
        ax.set_zlabel("MDS3", fontsize=8)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], c=cols, s=s)
        ax.set_ylabel("MDS2", fontsize=8)
    if annotate:
        for i, nm in enumerate(names):
            if dims >= 3:
                ax.text(coords[i, 0], coords[i, 1], coords[i, 2], nm, fontsize=label_size)
            else:
                ax.text(coords[i, 0], coords[i, 1], nm, fontsize=label_size)
    ax.set_xlabel("MDS1", fontsize=8)
    ax.tick_params(labelsize=7)
    return ax
