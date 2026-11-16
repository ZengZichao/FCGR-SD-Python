# Contributing

Thanks for your interest in FCGR-SD, a Python reproduction of Li et al. (2024),
*Exploring geometry of genome space via Grassmann manifolds* (The Innovation 5: 100677).

## Development setup

```bash
git clone https://github.com/ZengZichao/FCGR-SD-Python.git
cd FCGR-SD-Python
python -m venv .venv && source .venv/bin/activate
pip install -e ".[test]"          # numpy, pytest, scipy
pip install -e ".[plot,bench]"    # optional: matplotlib, biopython
```

## Checks before opening a pull request

```bash
python -m pytest tests -q
ruff check .
ruff format --check .
```

All tests must pass and ruff must report neither lint errors nor format
differences. Keep changes focused, add a test for every behaviour change, and
update both `docs/` and `docs/zh-CN/` when user-facing behaviour changes.

## Commit messages

Use short imperative summaries, for example `Add exact grid indexing test` or
`chore: update ruff configuration`.

## Reporting problems

Open an issue with a minimal reproducible example: a short FASTA excerpt, the
exact command you ran, and the observed versus expected output. For security
issues, follow [SECURITY.md](SECURITY.md) instead of opening a public issue.
