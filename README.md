# PdOx_CH4

## 1. Scripts

All automation tools for generating input and analyzing NEB results are provided in the scripts/ folder:

- generate_initial_final_states.py: Places CH4 molecule and generates CH3 + H final geometries with overlap checks.
- interpolate_neb_images.py: Generates NEB images between valid initial/final state pairs.
- run_neb_job.sh: Shell script template to submit NEB jobs to a cluster or local environment.

## 2. Optimized Surface Structures

The input_structures/ folder contains optimized surface models (in .vasp format) used to generate activation pathways:

- Pd.vasp, PdO.vasp, Pd-PdO.vasp, and PdOx.vasp

## 3. NEB Output Structures

The neb_outputs/ directory contains:

- Initial and final geometries for each surface site.
- Final NEB-converged transition states.
- Site-wise organization by surface type and site index.


## Citation

If you use this repository or scripts in your work, please cite:

> [**DOI Placeholder**] "Title of Manuscript", *Journal*, 2025.

---

For any questions or issues, feel free to open an issue or contact the authors.
