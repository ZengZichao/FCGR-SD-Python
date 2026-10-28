"""端到端流程编排：FCGR → 截断 SVD → 广义 Grassmann 距离 → NJ 树。"""

from __future__ import annotations

import math
import time
from collections import Counter
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from .distance import build_bases, build_svd_cache, distance_matrix_from_bases
from .grassmann import DEFAULT_EPS, DEFAULT_P
from .metrics import robinson_foulds
from .tree import Tree, neighbor_joining

__all__ = ["RunResult", "run_fcgr_sd", "suggest_n", "scan_truncation_p"]


def suggest_n(seq_length: int) -> int:
    """按论文经验给出 FCGR 分辨率参数：``n = ceil(log4(L))``。

    论文对 44 条 SARS-CoV-2 基因组（平均长 29,759）取 n = 8，
    对 30 条 HA 基因（平均长 1,718.8）取 n = 6，与此规则一致。
    """
    if seq_length < 4:
        return 1
    return max(1, int(math.ceil(math.log(seq_length, 4))))


@dataclass
class RunResult:
    """一次完整运行的结果。"""

    names: List[str]
    sequences: List[str]
    distances: np.ndarray
    bases: List[np.ndarray]
    ranks: np.ndarray
    tree: Optional[Tree]
    timings: Dict[str, float] = field(default_factory=dict)
    params: Dict[str, Any] = field(default_factory=dict)

    def summary(self) -> Dict[str, Any]:
        d = self.distances
        iu = np.triu_indices_from(d, k=1)
        return {
            "n_sequences": len(self.names),
            "n": self.params.get("n"),
            "p": self.params.get("p"),
            "alphabet": self.params.get("alphabet"),
            "truncation_rule": self.params.get("rule"),
            "mode": self.params.get("mode"),
            "normalize": self.params.get("normalize"),
            "eps": self.params.get("eps"),
            "method": self.params.get("method"),
            "clamp_negative": self.params.get("clamp_negative"),
            "rank_min": int(self.ranks.min()) if self.ranks.size else 0,
            "rank_median": float(np.median(self.ranks)) if self.ranks.size else 0.0,
            "rank_max": int(self.ranks.max()) if self.ranks.size else 0,
            "distance_min": float(d[iu].min()) if iu[0].size else 0.0,
            "distance_max": float(d[iu].max()) if iu[0].size else 0.0,
            "distance_mean": float(d[iu].mean()) if iu[0].size else 0.0,
            "timings": dict(self.timings),
        }


def run_fcgr_sd(
    names: Sequence[str],
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
    build_tree: bool = True,
    clamp_negative: bool = False,
    verbose: bool = False,
) -> RunResult:
    """执行完整的 FCGR-SD 流程。

    Parameters
    ----------
    names, sequences : sequence
        序列名与序列（顺序一致）。至少两条；名字必须唯一（作为叶标签）。
    n : int
        FCGR 分辨率参数，矩阵尺寸 ``2**n``。
    p : float
        截断阈值（奇异值累积贡献率）。
    alphabet : int
        碱基-顶点对应关系编号（1~6）。
    rule : {"ge", "le"}
        截断维数的判定方式，见 :func:`fcgrsd.grassmann.truncation_rank`。
    mode : {"exact", "float"}
        FCGR 网格定位方式，见 :func:`fcgrsd.fcgr.fcgr_matrix`。
    n_jobs : int
        成对距离计算的并行线程数。
    build_tree : bool
        是否用 NJ 建树。
    clamp_negative : bool
        是否把 NJ 的负枝长截断为 0（默认保留）。

    Returns
    -------
    RunResult
    """
    if len(names) != len(sequences):
        raise ValueError("names 与 sequences 长度不一致")
    if len(sequences) < 2:
        raise ValueError("至少需要两条序列（单条序列无法定义成对距离与系统发育树）")
    dup = [k for k, v in Counter(names).items() if v > 1]
    if dup:
        raise ValueError(f"序列名必须唯一，重复的有：{sorted(dup)[:5]}")

    timings: Dict[str, float] = {}
    t_all = time.perf_counter()

    if verbose:
        print("[Step 1] 计算 FCGR 与截断 SVD ...", flush=True)
    t0 = time.perf_counter()
    # 不经过 build_svd_cache：该缓存会为每条序列保留完整的 U（O(m·4^n)），
    # 在 n=8、10348 条时约 5.4 GB；build_bases 逐条处理，只保留截断基。
    bases = build_bases(
        sequences, n, p, alphabet=alphabet, rule=rule, mode=mode, normalize=normalize
    )
    timings["step1_2_fcgr_svd"] = time.perf_counter() - t0
    ranks = np.array([b.shape[1] for b in bases], dtype=np.int64)
    if verbose:
        print(
            f"         完成，用时 {timings['step1_2_fcgr_svd']:.3f}s；"
            f"截断维数 min/median/max = {ranks.min()}/{int(np.median(ranks))}/{ranks.max()}",
            flush=True,
        )

    if verbose:
        print("[Step 3] 计算成对广义 Grassmann 距离 ...", flush=True)
    t0 = time.perf_counter()
    D = distance_matrix_from_bases(bases, eps=eps, method=method, n_jobs=n_jobs)
    timings["step3_distance"] = time.perf_counter() - t0
    if verbose:
        m = len(sequences)
        print(
            f"         完成 {m * (m - 1) // 2} 对，用时 {timings['step3_distance']:.3f}s",
            flush=True,
        )

    tree: Optional[Tree] = None
    if build_tree:
        if verbose:
            print("[Step 4] Neighbor-Joining 建树 ...", flush=True)
        t0 = time.perf_counter()
        tree = neighbor_joining(names, D, clamp_negative=clamp_negative)
        timings["step4_nj"] = time.perf_counter() - t0
        if verbose:
            print(f"         完成，用时 {timings['step4_nj']:.3f}s", flush=True)

    timings["total"] = time.perf_counter() - t_all
    if verbose:
        print(f"[完成] 总用时 {timings['total']:.3f}s", flush=True)

    return RunResult(
        names=list(names),
        sequences=list(sequences),
        distances=D,
        bases=bases,
        ranks=ranks,
        tree=tree,
        timings=timings,
        params={
            "n": n,
            "p": p,
            "alphabet": alphabet,
            "rule": rule,
            "mode": mode,
            "normalize": normalize,
            "eps": eps,
            "method": method,
            "clamp_negative": clamp_negative,
        },
    )


def scan_truncation_p(
    names: Sequence[str],
    sequences: Sequence[str],
    n: int,
    *,
    ps: Optional[Sequence[float]] = None,
    reference_tree: Optional[Tree] = None,
    alphabet: int = 3,
    rule: str = "ge",
    mode: str = "exact",
    n_jobs: int = 1,
    verbose: bool = False,
) -> List[Dict[str, Any]]:
    """扫描截断参数 p（论文 §1.2.4 的调参流程）。

    论文对小数据集的做法是：在 ``0.01 ~ 0.99``（步长 0.01）上遍历 p，
    选使 NJ 树与 MAFFT 参考树的 RF 距离最小者。

    内存提示：为避免 99 次重复 SVD，本函数用 :func:`build_svd_cache` 为每条序列
    保留完整的左奇异矩阵 ``U``，占用 ``O(m·4^n)``；``n=8``、1 万条时约 5 GB。
    它面向论文规模的小数据集（30–44 条），大规模数据请直接用
    :func:`run_fcgr_sd`（逐条处理，只保留截断基）。

    Returns
    -------
    list of dict
        每项包含 ``p``、``rank_mean``、``rank_min``、``rank_max``，
        若给出 ``reference_tree`` 还包含 ``rf`` 与 ``rf_normalized``。
    """
    if ps is None:
        ps = np.round(np.arange(0.01, 1.0, 0.01), 2)
    cache = build_svd_cache(sequences, n, alphabet=alphabet, mode=mode)
    rows: List[Dict[str, Any]] = []
    for i, p in enumerate(ps):
        bases = build_bases(sequences, n, p, rule=rule, cache=cache)
        ranks = np.array([b.shape[1] for b in bases])
        D = distance_matrix_from_bases(bases, n_jobs=n_jobs)
        row: Dict[str, Any] = {
            "p": float(p),
            "rank_mean": float(ranks.mean()),
            "rank_min": int(ranks.min()),
            "rank_max": int(ranks.max()),
        }
        if reference_tree is not None:
            tree = neighbor_joining(names, D)
            row["rf"] = robinson_foulds(tree, reference_tree)
            row["rf_normalized"] = robinson_foulds(tree, reference_tree, normalize=True)
        rows.append(row)
        if verbose and (i + 1) % 10 == 0:
            print(f"  [scan-p] {i + 1}/{len(ps)}", flush=True)
    return rows
