"""Neighbor-joining, Newick I/O, and tree utilities.

The tree builder is the classical neighbor-joining algorithm of
Saitou & Nei (1987): at each step merge the pair minimizing the Q criterion,
update the distances, and stop when two nodes remain. The result is an
unrooted tree represented with a degree-2 root, written as standard Newick.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

import numpy as np

__all__ = ["Node", "Tree", "neighbor_joining", "read_newick"]


@dataclass
class Node:
    """树结点（叶子或内部结点）。"""

    name: Optional[str] = None
    length: Optional[float] = None
    children: List["Node"] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def leaves(self) -> List["Node"]:
        if self.is_leaf:
            return [self]
        out: List[Node] = []
        for c in self.children:
            out.extend(c.leaves())
        return out


def _fmt_length(x: float, digits: int) -> str:
    s = f"{x:.{digits}f}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    if s in ("", "-0", "-"):
        s = "0"
    return s


def _quote(name: str) -> str:
    if name is None or name == "":
        return ""
    unsafe = set(" \t\n()[],:;'\"")
    if any(ch in unsafe for ch in name):
        return "'" + name.replace("'", "''") + "'"
    return name


@dataclass
class Tree:
    """极简系统发育树容器。"""

    root: Node
    rooted: bool = False

    # -- 序列化 ------------------------------------------------------------
    def to_newick(self, digits: int = 10, include_lengths: bool = True) -> str:
        def rec(node: Node) -> str:
            if node.children:
                inner = ",".join(rec(c) for c in node.children)
                label = _quote(node.name or "")
                text = f"({inner}){label}"
            else:
                text = _quote(node.name or "")
            if include_lengths and node.length is not None:
                text += ":" + _fmt_length(node.length, digits)
            return text

        return rec(self.root) + ";"

    def __str__(self) -> str:  # pragma: no cover - 便于交互查看
        return self.to_newick()

    # -- 查询 --------------------------------------------------------------
    def leaf_names(self) -> List[str]:
        return [n.name or "" for n in self.root.leaves()]

    @property
    def n_leaves(self) -> int:
        return len(self.root.leaves())

    def splits(self) -> Set[FrozenSet[str]]:
        """非平凡二分划（unrooted bipartition）集合，用于 RF 距离。"""
        leaves = self.root.leaves()
        label = {}
        for i, node in enumerate(leaves):
            label[id(node)] = node.name if node.name else f"leaf{i}"
        total = len(leaves)
        all_names = frozenset(label.values())
        out: Set[FrozenSet[str]] = set()

        def rec(node: Node) -> Set[str]:
            if node.is_leaf:
                return {label[id(node)]}
            acc: Set[str] = set()
            for c in node.children:
                acc |= rec(c)
            if node is not self.root and 1 < len(acc) < total - 1:
                side = frozenset(acc)
                other = all_names - side
                key = side if (len(side), sorted(side)) <= (len(other), sorted(other)) else other
                out.add(key)
            return acc

        rec(self.root)
        return out

    def patristic_distances(self) -> Tuple[List[str], np.ndarray]:
        """返回叶子名与树上的路径距离矩阵（cophenetic distance）。"""
        leaves = self.root.leaves()
        names = [n.name if n.name else f"leaf{i}" for i, n in enumerate(leaves)]
        m = len(leaves)
        D = np.zeros((m, m))
        index = {id(n): i for i, n in enumerate(leaves)}

        def rec(node: Node) -> Dict[int, float]:
            """返回 {叶子全局下标: 到 node 的距离}。"""
            if node.is_leaf:
                return {index[id(node)]: 0.0}
            child_maps: List[Dict[int, float]] = []
            for c in node.children:
                w = c.length or 0.0
                child_maps.append({k: v + w for k, v in rec(c).items()})
            # 只有分属不同子结点的叶子对，其最近公共祖先才是当前结点
            for a in range(len(child_maps)):
                for b in range(a + 1, len(child_maps)):
                    for i, di in child_maps[a].items():
                        for j, dj in child_maps[b].items():
                            D[i, j] = D[j, i] = di + dj
            merged: Dict[int, float] = {}
            for cm in child_maps:
                merged.update(cm)
            return merged

        rec(self.root)
        return names, D

    def clusters(self) -> Set[FrozenSet[str]]:
        """有根树的内部分支簇（clade）集合。"""
        leaves = self.root.leaves()
        label = {}
        for i, node in enumerate(leaves):
            label[id(node)] = node.name if node.name else f"leaf{i}"
        total = len(leaves)
        out: Set[FrozenSet[str]] = set()

        def rec(node: Node) -> Set[str]:
            if node.is_leaf:
                return {label[id(node)]}
            acc: Set[str] = set()
            for c in node.children:
                acc |= rec(c)
            if node is not self.root and 1 < len(acc) < total:
                out.add(frozenset(acc))
            return acc

        rec(self.root)
        return out


def neighbor_joining(
    names: Sequence[str],
    distances: np.ndarray,
    *,
    clamp_negative: bool = False,
) -> Tree:
    """经典 Neighbor-Joining（Saitou & Nei 1987）。

    Parameters
    ----------
    names : sequence of str
        叶子名称，顺序须与距离矩阵一致。
    distances : ndarray
        形状 ``(m, m)`` 的对称距离矩阵。
    clamp_negative : bool
        是否把负枝长截断为 0。NJ 在距离不满足可加性时可能给出负枝长，
        NJ 在距离不满足可加性时可能给出负枝长；默认保留（``False``），
        需要非负枝长时设为 ``True``。

    Returns
    -------
    Tree
        无根树（以最后合并的结点为根）。
    """
    names = list(names)
    m = len(names)
    if m == 0:
        raise ValueError("至少需要一条序列")
    if len(set(names)) != len(names):
        raise ValueError("叶名必须唯一；重复的名字会让 Newick 与 RF 距离失去意义")
    D = np.array(distances, dtype=np.float64, copy=True)
    if D.shape != (m, m):
        raise ValueError(f"距离矩阵形状应为 ({m},{m})，收到 {D.shape}")
    if not np.all(np.isfinite(D)):
        raise ValueError("距离矩阵包含 NaN 或 Inf，无法建树")
    if not np.allclose(D, D.T, equal_nan=False):
        D = 0.5 * (D + D.T)
    np.fill_diagonal(D, 0.0)

    if m == 1:
        return Tree(root=Node(name=names[0]), rooted=False)

    nodes: List[Node] = [Node(name=nm) for nm in names]
    alive = np.ones(m, dtype=bool)
    size = m
    big = np.inf

    while size > 2:
        r = D.sum(axis=1)  # 已剔除的行/列被置零，不贡献
        Q = (size - 2) * D - r[:, None] - r[None, :]
        Q[~alive, :] = big
        Q[:, ~alive] = big
        np.fill_diagonal(Q, big)
        flat = int(np.argmin(Q))
        i, j = divmod(flat, m)
        if i > j:
            i, j = j, i

        dij = D[i, j]
        li = 0.5 * dij + (r[i] - r[j]) / (2.0 * (size - 2))
        lj = dij - li
        if clamp_negative:
            li = max(li, 0.0)
            lj = max(lj, 0.0)

        parent = Node(children=[nodes[i], nodes[j]])
        nodes[i].length = li
        nodes[j].length = lj
        nodes[i] = parent

        # 更新到其余活动结点的距离
        for k in range(m):
            if not alive[k] or k == i or k == j:
                continue
            dnew = 0.5 * (D[i, k] + D[j, k] - dij)
            D[i, k] = D[k, i] = dnew
        D[i, i] = 0.0

        alive[j] = False
        D[j, :] = 0.0
        D[:, j] = 0.0
        size -= 1

    rest = np.flatnonzero(alive)
    a, b = int(rest[0]), int(rest[1])
    d = float(D[a, b])
    half = d / 2.0
    if clamp_negative:
        half = max(half, 0.0)
    nodes[a].length = half
    nodes[b].length = half
    root = Node(children=[nodes[a], nodes[b]])
    return Tree(root=root, rooted=False)


_COMMENT_RE = None


def read_newick(text: str) -> Tree:
    """解析 Newick 字符串（支持引号名、枝长与 ``[...]`` 注释）。"""
    import re

    global _COMMENT_RE
    if _COMMENT_RE is None:
        _COMMENT_RE = re.compile(r"\[[^\]]*\]")
    s = _COMMENT_RE.sub("", text).strip()
    if s.endswith(";"):
        s = s[:-1]
    pos = 0
    n = len(s)

    def parse_node() -> Node:
        nonlocal pos
        node = Node()
        if pos < n and s[pos] == "(":
            pos += 1
            while True:
                node.children.append(parse_node())
                if pos < n and s[pos] == ",":
                    pos += 1
                    continue
                if pos < n and s[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"Newick 解析失败：位置 {pos} 附近 {s[max(0, pos - 20):pos + 20]!r}")
        # 标签
        if pos < n and s[pos] == "'":
            pos += 1
            buf = []
            while pos < n:
                if s[pos] == "'":
                    if pos + 1 < n and s[pos + 1] == "'":
                        buf.append("'")
                        pos += 2
                        continue
                    pos += 1
                    break
                buf.append(s[pos])
                pos += 1
            node.name = "".join(buf)
        else:
            start = pos
            while pos < n and s[pos] not in ",():;":
                pos += 1
            label = s[start:pos].strip()
            node.name = label or None
        # 枝长
        if pos < n and s[pos] == ":":
            pos += 1
            start = pos
            while pos < n and s[pos] not in ",():;":
                pos += 1
            token = s[start:pos].strip()
            if token:
                node.length = float(token)
        return node

    root = parse_node()
    return Tree(root=root, rooted=False)
