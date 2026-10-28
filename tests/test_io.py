"""FASTA 与距离矩阵读写的单元测试。"""

from __future__ import annotations

import numpy as np
import pytest

from fcgrsd.io import read_distance_matrix, read_fasta, write_distance_matrix


def test_distance_matrix_roundtrip(workdir) -> None:
    D = np.array([[0.0, 1.25, 2.5], [1.25, 0.0, 0.75], [2.5, 0.75, 0.0]])
    names = ["a", "b", "c"]
    path = workdir / "d.csv"
    write_distance_matrix(str(path), names, D)
    names2, D2 = read_distance_matrix(str(path))
    assert names2 == names
    assert np.allclose(D, D2, atol=1e-9)


def test_distance_matrix_header_mismatch_raises(workdir) -> None:
    """表头与首列标签不一致时必须报错，而不是静默采用表头。"""
    path = workdir / "bad.csv"
    path.write_text(",X,Y\nA,0,1\nB,1,0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="不一致"):
        read_distance_matrix(str(path))


def test_distance_matrix_non_square_raises(workdir) -> None:
    path = workdir / "bad2.csv"
    path.write_text(",a,b\nx,0,1\ny,1,0\nz,1,1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="方阵"):
        read_distance_matrix(str(path))


def test_read_fasta_rejects_empty_sequence(workdir) -> None:
    path = workdir / "empty.fasta"
    path.write_text(">a\n>b\nACGT\n", encoding="utf-8")
    with pytest.raises(ValueError, match="为空"):
        read_fasta(str(path))


def test_read_fasta_rejects_duplicate_names(workdir) -> None:
    path = workdir / "dup.fasta"
    path.write_text(">a\nACGT\n>a\nACGT\n", encoding="utf-8")
    with pytest.raises(ValueError, match="重复"):
        read_fasta(str(path))
