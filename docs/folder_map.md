# Where files belong

The active local project and GitHub use the same tracked folder structure.

| Folder | Contents | GitHub |
|---|---|---|
| `matlab/algorithms/` | RANSAC, PEARL, graph cuts, plane sampling and SVD | Yes |
| `matlab/config/` | MATLAB configuration function | Yes |
| `matlab/io/` | Dataset readers and grid/graph preparation | Yes |
| `matlab/evaluation/` | Metrics, matching, and output exports | Yes |
| `matlab/experiments/` | Experiment and tuning implementations | Yes |
| `matlab/scripts/` | Two experiment entry scripts and tuner entry script | Yes |
| `matlab/tests/` | Implementation checks | Yes |
| `matlab/visualization/` | Dataset and result-figure generators | Yes |
| `config/` | Portable confirmed parameter MAT file | Yes |
| `docs/` | Reproduction guide, validation notes, and citation | Yes |
| `data/` | Dataset setup instructions and archive checksum | Yes |
| `data/synpeb/`, `data/archives/` | Extracted dataset and downloaded archive | Local only |
| `results/reported/` | Compact CSV measurements used by the article | Yes |
| `results/runs/` | Full runs, per-case predictions, figures and checkpoints | Local only |
| `results/figures/` | Editable/generated MATLAB figure outputs | Local only |
| `paper/latex/` | Individual LaTeX sections, bibliography, IEEE class and figure PNGs | Yes |
| `paper/report/` | Reviewed report PDF | Yes |
| `paper/drafts/`, `paper/proposal/` | Earlier Word sections and proposal | Local only |
| `submission/` | Generated submission ZIPs and report copy | README only |
| `archive/` | Prior packages, old instructions, original code and parameter backup | Local only |

Run `addpath('matlab'); setup_project;` from the repository root before calling
functions. Entry scripts also initialize their paths when run directly. Saved
historical run configurations can contain old absolute paths; the current
figure generator resets the data root, and new runs use the new layout. To rerun
tuning after changing code or search settings, use a fresh tuning output folder.
