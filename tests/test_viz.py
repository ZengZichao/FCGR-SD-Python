"""可视化的单元测试（matplotlib 为可选依赖）。"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("matplotlib")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from fcgrsd.tree import neighbor_joining  # noqa: E402
from fcgrsd.viz import draw_tree, plot_fcgr  # noqa: E402


def test_plot_fcgr_matches_cgr_orientation() -> None:
    """FCGR 热图的横轴应为 CGR 第一坐标，故显示的是矩阵的转置。"""
    M = np.arange(16, dtype=float).reshape(4, 4)
    im = plot_fcgr(M, log=False)
    assert np.array_equal(np.asarray(im.get_array()), M.T)
    plt.close("all")


def test_draw_tree_zero_length() -> None:
    """全部枝长为 0 时不应产生退化的坐标轴。"""
    tree = neighbor_joining(["A", "B", "C"], np.zeros((3, 3)))
    ax = draw_tree(tree)
    lo, hi = ax.get_xlim()
    assert np.isfinite([lo, hi]).all()
    assert hi > lo
    plt.close("all")
