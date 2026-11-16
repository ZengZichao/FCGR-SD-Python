# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-10-28

### Added

- CGR/FCGR construction with closed-form integer grid indexing (`fcgrsd.cgr`,
  `fcgrsd.fcgr`).
- Truncated SVD bases, principal angles and the generalized Grassmann distance
  (`fcgrsd.grassmann`).
- Pairwise distance matrices in exact and float modes with optional threading
  (`fcgrsd.distance`).
- Neighbor-joining with Newick I/O (`fcgrsd.tree`); Robinson–Foulds distance and
  classical MDS (`fcgrsd.metrics`).
- End-to-end pipeline and the `fcgrsd` command-line interface (`run`, `scan-p`,
  `tree`, `mds`).
- 30 HA example dataset, precomputed p-scan results, figures, and English/Chinese
  documentation.
- 113 unit tests checked against independent references (SciPy, Biopython,
  brute-force enumeration, exact rational arithmetic).
- Continuous integration for Python 3.9–3.13, including lint and format checks.
