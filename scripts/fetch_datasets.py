"""Download sequences from NCBI by the accessions listed in a CSV table.

The repository does not redistribute most sequence data. This script takes your
own table (a CSV path, or a table name looked up in ``--data-dir``) and calls
NCBI E-utilities ``efetch`` in batches, writing FASTA (sequence name =
accession) plus a ``*_meta.csv`` with the detected group column.

The accession column is detected from headers containing ``accession``,
``gisaid``, ``header`` or ``id``; a ``Header`` value such as
``"LC710779.1 Avian coronavirus"`` is reduced to its first token.

Usage::

    # Explicit CSV path
    python scripts/fetch_datasets.py --table my_accessions.csv --out data/seqs.fasta

    # Look up TableS1.*ids.csv in a directory
    python scripts/fetch_datasets.py --table TableS1 --data-dir path/to/tables \
        --out data/ha.fasta

Note
----
GISAID accessions (e.g. ``EPI_ISL_*``) are usually not retrievable from NCBI;
they are recorded in ``*.missing.txt`` instead.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Dict, List, Sequence, Tuple

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


#: 登录号列的候选关键字（按优先级）。TableS2/5 用 GISAID ID，TableS3 的登录号在 Header 首 token。
_ACC_KEYS = ("accession", "gisaid", "header", "id")
#: 分组列的候选关键字
_GROUP_KEYS = ("gene", "type", "variant", "subgenus", "lineage", "pangolin", "clade")


def _pick_column(cols: Sequence[str], keys: Sequence[str]) -> str | None:
    for key in keys:
        for c in cols:
            if key in c.lower():
                return c
    return None


def _accession(value: str, column: str) -> str:
    v = value.strip()
    # TableS3 的 Header 形如 "LC710779.1 Avian coronavirus ..."，只取登录号
    if column.lower() == "header" or " " in v:
        return v.split()[0]
    return v


def load_table(path: str) -> Tuple[List[str], Dict[str, str]]:
    with open(path, "rt", encoding="utf-8-sig", errors="replace") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"{path} 为空")
    cols = list(rows[0].keys())
    acc_col = _pick_column(cols, _ACC_KEYS)
    if acc_col is None:
        raise ValueError(f"{path} 中找不到登录号列（候选关键字 {_ACC_KEYS}，实际列 {cols}）")
    group_col = _pick_column(cols, _GROUP_KEYS)
    accs: List[str] = []
    groups: Dict[str, str] = {}
    for r in rows:
        acc = _accession(r.get(acc_col, ""), acc_col)
        if not acc:
            continue
        accs.append(acc)
        groups[acc] = r.get(group_col, "").strip() if group_col else ""
    return accs, groups


def efetch_batch(accessions: Sequence[str], retries: int = 3) -> str:
    url = EFETCH + "?" + urllib.parse.urlencode(
        {
            "db": "nuccore",
            "id": ",".join(accessions),
            "rettype": "fasta",
            "retmode": "text",
        }
    )
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, TimeoutError) as exc:  # pragma: no cover
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"下载失败：{accessions[:2]}... ({last})")


def parse_fasta(text: str) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    name = None
    buf: List[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            header = line[1:].strip()
            name = header.split()[0]
            buf = []
        elif line.strip() and name is not None:
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Download sequences listed in an accession-table CSV")
    ap.add_argument("--table", required=True, help="CSV path, or a table name looked up in --data-dir")
    ap.add_argument(
        "--data-dir",
        default="examples/data",
        help="Directory searched when --table is a name (default: examples/data)",
    )
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch", type=int, default=20)
    ap.add_argument("--sleep", type=float, default=0.4)
    ap.add_argument("--max", type=int, default=0, help="Maximum number of records (0 = all)")
    args = ap.parse_args()

    table = args.table
    if not table.lower().endswith(".csv"):
        import glob

        pattern = os.path.join(args.data_dir, f"{args.table}.*ids.csv")
        hits = sorted(glob.glob(pattern))
        if not hits:
            print(
                f"No table matching {pattern!r}. Pass --table /path/to/table.csv "
                f"or put the CSV in --data-dir.",
                file=sys.stderr,
            )
            return 2
        table = hits[0]

    accs, groups = load_table(table)
    if args.max:
        accs = accs[: args.max]
    print(f"共 {len(accs)} 个登录号，来源 {table}")

    found: List[Tuple[str, str]] = []
    missing: List[str] = []
    for i in range(0, len(accs), args.batch):
        chunk = accs[i : i + args.batch]
        text = efetch_batch(chunk)
        recs = parse_fasta(text)
        got = {n for n, _ in recs}
        for acc in chunk:
            if acc in got:
                found.append((acc, dict(recs)[acc]))
            else:
                missing.append(acc)
        print(f"  {min(i + args.batch, len(accs))}/{len(accs)}  已获取 {len(found)}，缺失 {len(missing)}")
        time.sleep(args.sleep)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "wt", encoding="utf-8") as fh:
        for name, seq in found:
            fh.write(f">{name}\n")
            for j in range(0, len(seq), 70):
                fh.write(seq[j : j + 70] + "\n")
    meta = os.path.splitext(args.out)[0] + "_meta.csv"
    with open(meta, "wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "group"])
        for name, _ in found:
            w.writerow([name, groups.get(name, "")])
    if missing:
        miss_path = os.path.splitext(args.out)[0] + "_missing.txt"
        with open(miss_path, "wt", encoding="utf-8") as fh:
            fh.write("\n".join(missing) + "\n")
        print(f"缺失 {len(missing)} 条，已记录到 {miss_path}")

    print(f"已写出 {len(found)} 条序列：{args.out}")
    print(f"分组元数据：{meta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
