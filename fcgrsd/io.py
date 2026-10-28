"""输入输出：FASTA 读写与 NCBI/GISAID 登录号表解析。"""

from __future__ import annotations

import csv
import gzip
import os
from collections import Counter
from typing import Dict, List, Sequence, Tuple

__all__ = [
    "read_fasta",
    "write_fasta",
    "read_accession_table",
    "write_distance_matrix",
    "read_distance_matrix",
]


def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def read_fasta(path: str) -> Tuple[List[str], List[str]]:
    """读取 FASTA（支持 ``.gz`` 与多行序列）。

    Returns
    -------
    names : list of str
        序列名（``>`` 之后的第一个空白分隔字段）。
    sequences : list of str
        对应的序列字符串（保留原始大小写，未做过滤）。
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"找不到 FASTA 文件：{path}")
    names: List[str] = []
    seqs: List[str] = []
    buf: List[str] = []
    cur: str | None = None
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    names.append(cur)
                    seqs.append("".join(buf))
                header = line[1:].strip()
                cur = header.split()[0] if header else f"seq{len(names) + 1}"
                buf = []
            else:
                if cur is None:
                    raise ValueError(f"{path} 不是合法的 FASTA：首行不是 '>' 开头")
                buf.append(line.strip())
    if cur is not None:
        names.append(cur)
        seqs.append("".join(buf))
    if not names:
        raise ValueError(f"{path} 中没有解析到任何序列")
    empty = [n for n, s in zip(names, seqs) if not s]
    if empty:
        raise ValueError(f"以下序列为空：{empty[:5]}")
    if len(set(names)) != len(names):
        counts = Counter(names)
        dup = sorted(n for n, c in counts.items() if c > 1)[:5]
        raise ValueError(f"FASTA 中存在重复序列名 {dup}，请先去重或重命名")
    return names, seqs


def write_fasta(path: str, names: Sequence[str], sequences: Sequence[str], width: int = 70) -> None:
    """写出 FASTA 文件（每行 ``width`` 个碱基）。"""
    with open(path, "wt", encoding="utf-8") as fh:
        for name, seq in zip(names, sequences):
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


def read_accession_table(path: str) -> List[Dict[str, str]]:
    """读取带表头的登录号清单 CSV（每行至少含一个登录号列）。"""
    with open(path, "rt", encoding="utf-8-sig", errors="replace") as fh:
        return [dict(row) for row in csv.DictReader(fh)]


def write_distance_matrix(path: str, names: Sequence[str], D, digits: int = 10) -> None:
    """把距离矩阵写成带表头的 CSV。"""
    import numpy as np

    arr = np.asarray(D, dtype=np.float64)
    with open(path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + list(names))
        fmt = f"{{:.{digits}f}}"
        for name, row in zip(names, arr):
            writer.writerow([name] + [fmt.format(v) for v in row])


def read_distance_matrix(path: str) -> Tuple[List[str], "object"]:
    """读取 :func:`write_distance_matrix` 生成的 CSV。"""
    import numpy as np

    with open(path, "rt", encoding="utf-8-sig", newline="") as fh:
        rows = [r for r in csv.reader(fh) if r]
    if not rows:
        raise ValueError(f"{path} 为空")
    header = rows[0][1:]
    row_names = [r[0] for r in rows[1:]]
    mat = np.array([[float(x) for x in r[1:]] for r in rows[1:]], dtype=np.float64)
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        raise ValueError(f"距离矩阵必须是方阵，收到形状 {mat.shape}")
    if not np.all(np.isfinite(mat)):
        raise ValueError("距离矩阵包含 NaN 或 Inf")
    if header and len(header) == mat.shape[1]:
        # 表头与首列标签都应存在；二者不一致时以报错代替静默覆盖
        if any(row_names) and row_names != header:
            raise ValueError("距离矩阵的表头与首列标签不一致，无法确定叶名")
        names = header
    else:
        names = row_names
    if len(names) != mat.shape[0]:
        raise ValueError(f"叶名数量 {len(names)} 与矩阵维数 {mat.shape[0]} 不一致")
    return names, mat
