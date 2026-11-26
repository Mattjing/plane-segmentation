# Validation of the organized source package

Verified in MATLAB R2026a on 2026-11-25 after restructuring:

- All 21 MATLAB files (20 original sources plus `setup_project.m`) passed the parse-error check.
- The setup function discovers all folders from the repository root.
- Full-mode configuration resolves `data/synpeb/`, `config/selected_parameters.mat`, and `results/runs/` correctly, with the original frozen settings.
- `test_synpeb_implementation` passed 25 exhaustive binary expansion checks, SVD, reference matching, and PEARL energy-descent checks.
- One real test scan loaded 3,969 points and a geometry-based graph with 7,732 edges from the new data path. RANSAC followed by one outer iteration and one expansion sweep of hybrid refinement completed with non-increasing energy.
- The article compiled to 17 pages with 14 references, no unresolved citations/references, and no overfull boxes.

No full experiment or hyperparameter search was repeated during reorganization.
The recorded CSV measurements and segmentation algorithms are unchanged.
The distributed parameter MAT file still contains only `selected`; original
search state and historical local paths remain in the local archive.

Version `v1.0.0` preserves the earlier flat layout. Version `v1.1.0` introduces
the organized layout and renamed public URL without rewriting that history.
