# Energy-Based Multi-Plane Segmentation of 3D Point Clouds

MATLAB source code and reproduction instructions for the ECE 503 course project by **Matt (Haocun) Jing**, University of Victoria. This repository accompanies a comparison experiment; it does not propose a new segmentation algorithm or reproduce every component of the original PEARL software.

**Public repository:** <https://github.com/Mattjing/ece503-plane-segmentation>  
**Report code version:** [`v1.0.0`](https://github.com/Mattjing/ece503-plane-segmentation/tree/v1.0.0)

## What is compared?

| Pipeline | Initialization and optimization |
|---|---|
| Sequential RANSAC | Repeatedly sample three points, retain a plane with strong consensus, refit by SVD, assign and remove its support, and continue on the remaining points. |
| Sampled-proposal PEARL | Generate and screen a bank of candidate planes, initialize labels by unary fitting costs, then alternate alpha-expansion assignment updates and SVD plane refitting. |
| RANSAC-initialized PEARL (hybrid) | Transfer the **exact planes and assignments** from the standalone RANSAC result, then refine them with the same PEARL optimizer. No additional candidate planes are introduced. |

A separate ablation compares sampled PEARL **with versus without spatial regularization**, using the same initial plane bank. This is not a tuned-versus-default ablation.

There is no supervised weight-training stage. Hyperparameters are selected on annotated training scans and frozen before testing; plane parameters and point assignments are estimated independently for every input scan. Reference annotations are used for parameter selection and evaluation, not as fitting features.

## Requirements

- **MATLAB R2026a** was used for the completed experiments and implementation checks. Earlier versions have not been verified.
- The implementation uses MATLAB functions including `digraph`, `maxflow`, `matchpairs`, `readtable`, `tiledlayout`, and `exportgraphics`.
- No external graph-cut MEX library, Computer Vision Toolbox, or Statistics and Machine Learning Toolbox is required by this implementation.
- Download SynPEB separately. The dataset and its archive are **not redistributed** in this repository. Consult the [dataset authors' PPE repository](https://github.com/acschaefer/ppe) and retain the dataset's accompanying license and attribution.

## Installation and dataset layout

```bash
git clone https://github.com/Mattjing/ece503-plane-segmentation.git
cd ece503-plane-segmentation
git checkout v1.0.0
```

Extract the downloaded dataset so that `synpeb` is directly beside the MATLAB scripts:

```text
ece503-plane-segmentation/
  run_synpeb_project.m
  run_initialization_comparison.m
  selected_parameters.mat
  synpeb/
    train/
      gt/pc_01.pcd ... pc_10.pcd
      var1/ ... var4/
    test/
      gt/pc_01.pcd ... pc_30.pcd
      var1/ ... var4/
```

The actual split layout has separate `gt`, `var1`, `var2`, `var3`, and `var4` directories. Avoid an extra nested `synpeb/synpeb` directory. The default data location is set in `synpeb_config.m`; change `cfg.dataRoot` there if necessary.

The source PCD files contain organized 500-by-500 scans. Their `intensity` field encodes reference IDs rather than an input reflectance feature. Each measurement variant is evaluated using **its own** intensity labels; labels are not blindly copied from a corresponding `gt` file. IDs below 10 are excluded from reference-plane metrics, while valid coordinates remain available to the fitting methods. At stride eight, a complete scan retains 63-by-63 positions before invalid-coordinate filtering.

## Quick implementation check

Set MATLAB's Current Folder to this repository and run:

```matlab
addpath(pwd);
test_synpeb_implementation;
```

This check does not require the downloaded dataset. It compares 25 small alpha-expansion moves against exhaustive binary minimization and checks SVD fitting, label matching, rejection and fragmentation behavior, and PEARL energy descent. It verifies the implemented binary moves; it does not establish global optimality of multi-model fitting.

For a small dataset run, edit `mode = 'full'` to `mode = 'demo'` in `run_synpeb_project.m`, then run:

```matlab
run_synpeb_project
```

Demo mode uses training scan 1, variants `var1` and `var4`, seed 7, and untuned demonstration defaults. **Demo output does not reproduce the report's final test scores.** Restore `mode = 'full'` for the following experiments.

## Reproduce the final experiments

`selected_parameters.mat` contains the confirmed **selected** configuration and its selection metadata. It is exported from the completed local tuning result with the large internal search state removed. It does not contain a learned classifier or reusable scene-specific planes. The included settings allow running the final experiments without repeating the lengthy search.

### 1. Three-way comparison

Keep `mode = 'full'` in `run_initialization_comparison.m` and run:

```matlab
run_initialization_comparison
```

This evaluates sequential RANSAC, sampled PEARL, and the hybrid on 30 test scans, all four measurement variants, and seeds `[7 19 43]`: **360 records per method**. The hybrid reuses the baseline RANSAC solution rather than rerunning initialization. Its reported total time includes the RANSAC initialization cost once plus refinement.

Outputs are written to:

```text
results/initialization_comparison/full/<timestamp>/
```

The folder includes `initialization_metrics.csv`, `initialization_accuracy_summary.txt`, `paired_refinement_changes.csv`, saved configuration and checkpoints, and per-case predictions, metrics, confusion matrices, accuracy text, segmentations, and PEARL energy histories.

### 2. Spatial-regularization ablation

Keep `mode = 'full'` in `run_synpeb_project.m` and run:

```matlab
run_synpeb_project
```

This runner evaluates RANSAC, sampled PEARL, and sampled PEARL with `lambdaRatio = 0`. The no-smoothness condition reuses the sampled proposal bank. Results are saved under `results/full/<timestamp>/`, including `metrics_per_case.csv` and `accuracy_summary.txt`.

The ablation uses the selected configuration with only spatial coupling removed; its unregularized condition is **not separately retuned**. Label-visit random trajectories are not explicitly locked between conditions. Compare runtimes within their own run, because the three-way and ablation timing sessions differed.

### 3. Generate dataset and results illustrations

Dataset illustrations use fixed training examples and the included confirmed settings:

```matlab
generate_dataset_figures;
```

After completing both experiments, supply the actual output paths from your runs:

```matlab
projectDir = pwd;
comparisonDir = fullfile(projectDir, 'results', ...
    'initialization_comparison', 'full', '<comparison-timestamp>');
ablationDir = fullfile(projectDir, 'results', 'full', '<ablation-timestamp>');
outputDir = fullfile(projectDir, 'results', 'report_figures');
generate_results_figures(projectDir, outputDir, comparisonDir, ablationDir);
```

Replace the timestamp placeholders with real folder names. `generate_results_figures` reads saved predictions and energy histories as well as CSV metrics; the compact CSVs included in `reported_results/` alone are insufficient to reconstruct the segmentation panels. Full generated results are not committed. The function's two-argument form retains the historical run-directory defaults used for the report. The optional directory arguments are a portability improvement; the segmentation algorithms are unchanged.

## Rerun hyperparameter selection from scratch

The included confirmed settings are sufficient for evaluation. To independently repeat their selection:

1. In `tune_synpeb_project.m`, set `phase = 'screen'` and choose a fresh `options.outputRoot` if algorithm code or search settings have changed. Run `tune_synpeb_project`.
2. After screening completes, set `phase = 'confirm'` and run it again with the same output root. The confirmation stage requires the screening shortlist.
3. Only completed confirmation writes the final `selected_parameters.mat`. Back up the included file before intentionally replacing it.
4. Run the full experiments with the resulting fixed settings. Do not select settings from held-out test performance.

Screening uses training scans 1–2, variants `var1` and `var4`, seed 7. The staged search varies proposal count, plane cap, local radius, 25 threshold–smoothness combinations, and minimum support. The three highest-ranked configurations are confirmed on disjoint training scans 3–10, all variants, and three seeds: 96 records per configuration and method. Selection maximizes the equally weighted mean matched IoU of RANSAC and PEARL, breaking ties by mean F-score and then summed runtime. This is a heuristic staged search, not an exhaustive search over every parameter combination. Checkpoints reuse completed configurations.

## Implemented objective and iteration budgets

For point labels $L_i$ and unit-normal planes $\theta_k=(n_k,d_k)$, the implemented basic PEARL energy is

$$E(L,\Theta)=\sum_i D_i(L_i;\Theta)+\lambda\sum_{(i,j)\in\mathcal E}w_{ij}[L_i\ne L_j].$$

For a plane label $k>0$, $D_i(k;\Theta)=(n_k^\top x_i+d_k)^2$; the outlier label 0 has cost $D_i(0;\Theta)=\tau^2$. The implementation sets $\lambda=\eta\tau^2$ and optimizes the equivalent normalized energy $E/\tau^2$. **There is no model-count penalty.**

Alpha-expansion updates assignments with planes fixed. Each binary move is solved by `maxflow`; full multi-label optimization remains approximate. SVD refitting updates each supported plane with assignments fixed and reduces the geometric term, while leaving the spatial term unchanged. Unused plane labels are removed only when producing the final output. The hybrid cannot create more plane labels than its initial RANSAC bank.

| Setting | Report configuration |
|---|---:|
| Organized-grid stride | 8 |
| Distance threshold / outlier scale, `distanceThreshold` | 0.04 |
| Smoothness ratio, `lambdaRatio` | 0.4 |
| Physical smoothness weight, $\lambda$ | 0.00064 |
| Outlier cost, $\tau^2$ | 0.0016 |
| Sampled PEARL proposal attempts | 1,000 |
| Maximum plane bank | 150 |
| Minimum proposal/extraction support | 6 |
| RANSAC attempts per extraction | 150 |
| Local sampling fraction / radius | 0.7 / 8 sampled-grid pixels |
| Candidate duplicate angle / offset limits | 3 degrees / 0.015 |
| Graph distance scale / maximum edge length | 0.15 / 0.30 |
| PEARL outer iterations / expansion sweeps | 15 / 3 |
| Normalized energy tolerance | $10^{-9}$ |

An expansion **sweep visits every plane label and the outlier label**, so three sweeps do not mean three individual graph cuts. With $M$ planes, the configured upper bound is $15\times3\times(M+1)$ binary moves; early stopping can reduce this. RANSAC permits up to three SVD refit-and-reclassification passes per accepted candidate.

## Accuracy, overlap, and saved outputs

Plane IDs are arbitrary within each case, so reference and predicted IDs are matched **one-to-one by maximum point overlap** before calculating point accuracy. Accuracy is the fraction of evaluable reference-plane points assigned to the matched prediction. Rejected points and unmatched predicted labels count as errors. Mean matched IoU uses that same matching, assigning zero to unmatched reference planes. Instance precision, recall, and F-score use a separate one-to-one matching with eligible IoU at least 0.5.

Each case saves a count-based and row-normalized confusion-matrix figure (`*_confusion.png` and `.fig`), numerical accuracy and other metrics (`*_accuracy.txt`), label and confusion CSVs, a MAT result, and segmentation figures. PEARL additionally saves an energy history. Complete label order is retained for dense confusion matrices. No aggregate confusion matrix is pooled across scans, because their plane IDs have different meanings.

Algorithm timings exclude loading, graph construction, evaluation matching, and figure export. Repeated seeds and variants of the same scan are correlated; reported standard deviations describe variation across records and are not confidence intervals from 360 independent scenes.

## Recorded report results

The CSVs in `reported_results/` are copies of the completed experiments, not newly generated measurements. Accuracy below is a record mean; times come only from the three-way comparison session.

| Method | Matched point accuracy (%) | Mean matched IoU | Mean F-score | Mean total time (s) |
|---|---:|---:|---:|---:|
| Sequential RANSAC | 69.89 | 0.2086 | 0.2257 | 0.376 |
| Sampled PEARL | 84.26 | 0.4554 | 0.4003 | 19.620 |
| Hybrid | 77.98 | 0.2836 | 0.3443 | 2.512 |

The hybrid improves point accuracy in 357 of 360 paired records and worsens it in three. It is about 7.81 times faster than sampled PEARL in that run, while sampled PEARL achieves higher accuracy and overlap. Removing spatial regularization reduces sampled PEARL mean accuracy from 84.26% to 11.27% under the selected fixed configuration. Energy descent does not guarantee improved reference-label accuracy.

Saved source-run identifiers are `20261124_022344_155` (three-way comparison) and `20261124_062337_283` (spatial ablation). Runtime depends on hardware and session conditions; a fresh run need not reproduce the exact recorded times. Reproducing the reference scores also requires the same dataset, sampling, seeds, settings, and software behavior.

## Source map

| Files | Purpose |
|---|---|
| `run_synpeb_project.m`, `run_synpeb_experiment.m` | Main RANSAC/PEARL evaluation and spatial-prior ablation |
| `run_initialization_comparison.m`, `compare_pearl_initializations.m` | Three-way initialization comparison and paired refinement |
| `sequential_ransac.m`, `sample_plane.m`, `fit_plane_svd.m` | Consensus extraction, triplet proposals, orthogonal least squares |
| `pearl_segment.m`, `pearl_refine_ransac.m`, `alpha_expansion_move.m` | Sampled and hybrid PEARL, binary graph cuts |
| `synpeb_config.m`, `tune_synpeb_project.m`, `synpeb_tuning_search.m` | Configuration and staged training-split parameter selection |
| `read_synpeb_pcd.m`, `load_synpeb_case.m` | PCD/LZF decoding, grid sampling, neighborhood construction |
| `evaluate_planes.m`, `save_case_results.m` | Matching, metrics, confusion matrices, accuracy text, figures |
| `generate_dataset_figures.m`, `generate_results_figures.m` | MATLAB illustrations used in the report |
| `test_synpeb_implementation.m` | Small exhaustive-cut, geometry, matching, and energy checks |
| `selected_parameters.mat` | Frozen confirmed configuration; only the `selected` variable is distributed |

## Scope and attribution

This is a course implementation of the basic PEARL pipeline, not the original authors' implementation. The comparison does not equalize proposal-bank size or computational budgets, and the hybrid is not separately tuned. SynPEB is synthetic; the results do not establish general robustness on real scenes or arbitrary unorganized point clouds. Stanford S3DIS semantic/object annotations are not directly equivalent to explicit plane-instance labels.

Dataset ownership and its accompanying terms remain with its authors. This repository excludes the SynPEB files, downloaded third-party source, local filesystem paths, old script backups, and generated experiment directories. Public availability of this code is separate from any redistribution rights to the dataset.

Method references:

1. M. A. Fischler and R. C. Bolles, “Random sample consensus: A paradigm for model fitting with applications to image analysis and automated cartography,” *Communications of the ACM*, 24(6), 381–395, 1981. [DOI](https://doi.org/10.1145/358669.358692).
2. H. Isack and Y. Boykov, “Energy-based geometric multi-model fitting,” *International Journal of Computer Vision*, 97(2), 123–147, 2012. [DOI](https://doi.org/10.1007/s11263-011-0474-7).
3. Y. Boykov, O. Veksler, and R. Zabih, “Fast approximate energy minimization via graph cuts,” *IEEE TPAMI*, 23(11), 1222–1239, 2001. [DOI](https://doi.org/10.1109/34.969114).
4. A. Schaefer, J. Vertens, D. Büscher, and W. Burgard, “A maximum likelihood approach to extract finite planes from 3-D laser scans,” *ICRA*, 72–78, 2019. [DOI](https://doi.org/10.1109/ICRA.2019.8794318).

## Cite the code

GitHub also provides citation metadata through [`CITATION.cff`](CITATION.cff). Cite the repository README and code version using:

```bibtex
@misc{jing2026planes,
  author       = {Jing, Matt (Haocun)},
  title        = {{Energy-Based Multi-Plane Segmentation of 3D Point Clouds: MATLAB Comparison Code}},
  year         = {2026},
  howpublished = {GitHub repository, version v1.0.0},
  url          = {https://github.com/Mattjing/ece503-plane-segmentation},
  note         = {Source code and reproduction README. Accessed: 2026-11-25}
}
```

The code submission accompanies the article separately because full listings are lengthy. Check the course submission requirements when packaging the final report and source archive.
