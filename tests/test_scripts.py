"""脚本层（数据表列识别）的单元测试。"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from fetch_datasets import load_table  # noqa: E402


@pytest.mark.parametrize(
    "header, value, group_header, group_value",
    [
        ("NCBI Accession Number", "KM409247.1", "Gene type", "H1N1"),
        ("GISAID ID", "EPI_ISL_10113812", "Variants", "Beta"),
        ("Header", "LC710779.1 Avian coronavirus, complete genome", "Subgenus", "Igacovirus"),
        ("Accession ID", "EPI_ISL_14466156", "Type", "Omicron BA.5"),
        ("GISAID_ID", "EPI_ISL_14593879", "Pangolin", "B.1.617.2"),
    ],
)
def test_load_table_column_variants(workdir, header, value, group_header, group_value) -> None:
    import csv

    path = workdir / "t.csv"
    with open(path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([header, group_header])
        w.writerow([value, group_value])
    accs, groups = load_table(str(path))
    expected = value.split()[0]
    assert accs == [expected]
    assert groups[expected] == group_value


def test_load_table_without_accession_column_raises(workdir) -> None:
    path = workdir / "bad.csv"
    path.write_text("Foo,Bar\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="登录号"):
        load_table(str(path))
