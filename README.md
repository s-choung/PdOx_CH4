# PdOx_CH4

## 1. Scripts

All automation tools for generating input and analyzing NEB results are provided in the `scripts/` folder:

- `generate_initial_final_states.py`: Places CH₄ molecule and generates CH₃ + H final geometries with overlap checks.
- `interpolate_neb_images.py`: Generates NEB images between valid initial/final state pairs.
- `run_neb_job.sh`: Shell script template to submit NEB jobs to a cluster or local environment.
- `visualize_trajectories.py`: Generates overlayed CH₄ activation geometries for each model.

## 2. Optimized Surface Structures

The `input_structures/` folder contains optimized surface models (in `.cif` format) used to generate activation pathways:

- `Pd.cif`, `PdO.cif`, `Pd-PdO.cif`, and `PdOx.cif`

## 3. NEB Output Structures

The `neb_outputs/` directory contains:

- Initial and final geometries for each surface site.
- Interpolated NEB images.
- Final NEB-converged transition states and output logs.
- Site-wise organization by surface type and site index.

## 4. Visualization

Overlay plots of all CH₄ activation trajectories are available in the `visualizations/` folder:

- These `.xyz` files can be opened in visualization tools (e.g., VESTA, ASE GUI, OVITO) to examine the distribution of reaction geometries across all surface sites.
- Color code: Pd (turquoise), O (red), Ce (beige), C (gray), H (white).

## Citation

If you use this repository or scripts in your work, please cite:

> [**DOI Placeholder**] "Title of Manuscript", *Journal*, 2025.

---

For any questions or issues, feel free to open an issue or contact the authors.
