# Recorded measurements

These files preserve the completed report's measurements without committing the dataset or large prediction/figure directories.

- `three_way_metrics.csv`: `initialization_metrics.csv` from comparison run `20261124_022344_155`; 360 records per method, three methods.
- `spatial_ablation_metrics.csv`: `metrics_per_case.csv` from run `20261124_062337_283`; RANSAC, PEARL, and PEARL without spatial regularization.
- `screen_results.csv`, `screen_shortlist.csv`, `confirm_results.csv`: completed training-split search results and shortlist.

The two experiment sessions have different runtime measurements. Compare runtimes within each session. Scan variants and repeated random seeds are correlated records, not independent scenes. Full segmentation visualization requires the saved per-case predictions from a reproduced run, in addition to the dataset and these CSVs.
