# Plane Segmentation

MATLAB comparison of sequential RANSAC, sampled-proposal PEARL, and
RANSAC-initialized PEARL on SynPEB. The repository accompanies the course report
*Energy-Based Multi-Plane Segmentation of 3D Point Clouds* by Matt (Haocun) Jing.

**Repository:** https://github.com/Mattjing/plane-segmentation  
**Current organized code version:** `v1.1.0` (same segmentation algorithms and reported measurements as `v1.0.0`).

## Folder guide

```text
matlab/
  setup_project.m      add all MATLAB folders to the path
  algorithms/          plane sampling, SVD, RANSAC, PEARL, alpha-expansion
  config/              configuration function
  io/                  PCD readers and point/graph preparation
  evaluation/          matching, metrics, figures and accuracy exports
  experiments/         experiment and tuning implementations
  scripts/             experiment and tuning entry points
  tests/               numerical implementation checks
  visualization/       dataset and result illustrations
config/                confirmed parameters
data/                  dataset instructions and checksum; data downloaded separately
docs/                  reproduction guide, validation, folder map, citation
results/reported/      compact recorded measurements
paper/latex/           editable article, sections, references, figures and IEEE class
paper/report/          compiled report PDF
submission/            instructions for local submission packages
```

Each source file is stored individually in its designated folder. The
[complete folder map](docs/folder_map.md) also documents local-only data, full
experiment outputs, earlier drafts, and archives.

## Quick start

Tested with MATLAB R2026a. Download SynPEB separately following
[data/README.md](data/README.md), then open this repository root in MATLAB:

```matlab
addpath('matlab');
setup_project;
test_synpeb_implementation;    % does not need the dataset
run_initialization_comparison % full three-way comparison
run_synpeb_project            % full spatial-prior ablation
```

The two experiment scripts default to full evaluation with the confirmed
settings in [config/selected_parameters.mat](config/selected_parameters.mat).
For a small training-only demonstration, change `mode` to `'demo'` in the
relevant entry script. Demo defaults do not reproduce the report's scores.

- [Full reproduction instructions and results](docs/reproducing.md)
- [Validation checks](docs/validation.md)
- [Individual paper source](paper/latex/main.tex)
- [Compiled report](paper/report/plane_segmentation_report.pdf)
- [Recorded metrics](results/reported/)
- [BibTeX code citation](docs/citation.bib) and [GitHub citation metadata](CITATION.cff)

This is a comparison implementation of existing methods. Hyperparameters are
selected on training scans; there are no classifier weights learned across
scans. The full dataset, old drafts, and large generated runs are kept outside
Git. Dataset ownership and accompanying terms remain with its authors.
