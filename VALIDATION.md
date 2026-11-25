# Validation of the public source package

Verified in MATLAB R2026a on 2026-11-25 before publication:

- All 20 MATLAB source files passed the parse-error check.
- `test_synpeb_implementation` passed 25 alpha-expansion comparisons against exhaustive binary minimization, SVD fitting, reference matching, and PEARL energy-descent checks.
- Full-mode configuration loaded the distributed parameter file and recovered threshold 0.04, smoothness ratio 0.4, 1,000 proposals, at most 150 planes, 15 outer iterations, three sweeps, and seeds `[7 19 43]`.
- The distributed MAT file contains only the `selected` variable. The local search checkpoint state is not included.

The source package retains the completed experiment's algorithms. The results-figure helper additionally accepts optional comparison and ablation directory arguments for portability. Full dataset experiments were not rerun as part of packaging; the provided metric CSVs are copied from the completed runs identified in the README.
