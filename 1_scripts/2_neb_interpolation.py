"""
Generate initial and final geometries for CH₄ activation NEB runs on
Pd–CeO₂ slabs with varying oxygen content.

The script:
1. Loads the best‐energy slab models (Pd45_O{0,20,45}) produced during GA
   optimisation.
2. Identifies surface Pd atoms.
3. Builds an initial state with a CH₄ molecule placed ~3 Å above each
   surface Pd.
4. Builds a final state where CH₃ is adsorbed on the same Pd and the
   dissociated H is attached to the nearest O (preferred) or Pd.
5. Performs relaxed geometric sanity checks to avoid unphysical overlaps.
6. Saves `*_initial.vasp` and `*_final.vasp` files ready for NEB in an
   organised directory tree.

All heavy external dependencies (ASE, NumPy, SciPy) are already present
in typical atomistic Python environments.

Author: Seokhyun Choung
Licence: MIT
"""

import os
import re
import numpy as np
from ase.io import read, write
from ase.neb import NEB
from ase.optimize import LBFGS
from ase.io.vasp import read_vasp, write_vasp
from ase import Atoms
from ase.neighborlist import neighbor_list
from scipy.spatial.transform import Rotation as R

# CONFIGURATION
BASE_DIR = './neb_input_relaxed'
OUTPUT_DIR = './neb_interpolated'
os.makedirs(OUTPUT_DIR, exist_ok=True)

O_NUMS = [0, 20, 45]  # slab oxidation states
CLEARANCE_CH4 = 3.0  # Å C–Pd distance in initial state
CLEARANCE_CH3 = 2.0  # Å C–Pd distance in final state
SURFACE_THRESHOLD_Z = 6.0  # Å Pd considered “surface Pd”

# Distance cut-offs for overlap test (Å)
H_MIN_DIST = 0.8
DEFAULT_MIN_DIST = 1.8


def detect_ch3_and_metal_bond(atoms):
    """
    Detect CH₃ and the metal atom it's bonded to in the final structure

    Returns:
    - ch3_indices: Indices of C and H atoms in CH3
    - metal_idx: Index of the metal atom bonded to C
    - dissociated_h_idx: Index of the dissociated H atom
    """
    symbols = atoms.get_chemical_symbols()
    c_indices = [i for i, sym in enumerate(symbols) if sym == 'C']
    
    # Find CH3 group (carbon with 3 hydrogens)
    ch3_indices = []
    metal_idx = None
    dissociated_h_idx = None
    
    for c_idx in c_indices:
        c_pos = atoms.positions[c_idx]
        
        # Find hydrogens close to this carbon (CH bond length is typically ~1.1 Å)
        h_indices = [i for i, sym in enumerate(symbols) if sym == 'H' and 
                     np.linalg.norm(atoms.positions[i] - c_pos) < 1.5]
        
        if len(h_indices) == 3:  # If this C has 3 H atoms, it's likely CH3
            ch3_indices = [c_idx] + h_indices
            
            # Find metal atom bonded to carbon (looking for Pd within bonding distance)
            metal_indices = [i for i, sym in enumerate(symbols) if sym == 'Pd' and 
                            np.linalg.norm(atoms.positions[i] - c_pos) < 2.5]
            
            if metal_indices:
                distances = [np.linalg.norm(atoms.positions[i] - c_pos) for i in metal_indices]
                metal_idx = metal_indices[np.argmin(distances)]
            
            # Find the dissociated H (should be on the surface and away from C)
            all_h_indices = [i for i, sym in enumerate(symbols) if sym == 'H']
            remaining_h = [h for h in all_h_indices if h not in h_indices]
            
            if remaining_h:
                distances = [np.linalg.norm(atoms.positions[h] - c_pos) for h in remaining_h]
                dissociated_h_idx = remaining_h[np.argmax(distances)]
            
            break
    
    return ch3_indices, metal_idx, dissociated_h_idx


def get_site_numbers(folder):
    """Get site numbers from the folder"""
    files = os.listdir(folder)
    site_numbers = set()
    for fname in files:
        match = re.search(r'_site(\d+)_initial.vasp', fname)
        if match:
            site_numbers.add(int(match.group(1)))
    return sorted(site_numbers)


def generate_initial_ch4_structure(final_structure):
    """
    Generate initial CH4 structure from the final structure by back-propagation.
    
    Returns:
    - initial_structure: Generated initial structure with CH4 molecule
    - ch4_indices: Indices of atoms in CH4 [C, H1, H2, H3, H4]
    """
    initial = final_structure.copy()
    
    ch3_indices, metal_idx, diss_h_idx = detect_ch3_and_metal_bond(final_structure)
    
    if not ch3_indices or metal_idx is None or diss_h_idx is None:
        raise ValueError("Could not detect CH3-Pd bond or dissociated H in final structure")
    
    c_idx = ch3_indices[0]
    ch3_h_indices = ch3_indices[1:]
    
    c_pos = final_structure.positions[c_idx]
    metal_pos = final_structure.positions[metal_idx]
    diss_h_pos = final_structure.positions[diss_h_idx]
    ch3_h_pos = final_structure.positions[ch3_h_indices]
    
    c_pd_vector = metal_pos - c_pos
    c_pd_distance = np.linalg.norm(c_pd_vector)
    c_pd_unit_vector = c_pd_vector / c_pd_distance
    
    new_c_pos = metal_pos - c_pd_unit_vector * (c_pd_distance + 1.5)
    initial.positions[c_idx] = new_c_pos
    
    ch3_vectors = ch3_h_pos - c_pos
    ch3_distances = np.linalg.norm(ch3_vectors, axis=1)
    
    for i, h_idx in enumerate(ch3_h_indices):
        initial.positions[h_idx] = new_c_pos + ch3_vectors[i]
    
    avg_vector = np.mean(ch3_vectors, axis=0)
    avg_distance = np.mean(ch3_distances)
    opposite_vector = -avg_vector
    opposite_unit_vector = opposite_vector / np.linalg.norm(opposite_vector)
    new_h_pos = new_c_pos + opposite_unit_vector * avg_distance
    initial.positions[diss_h_idx] = new_h_pos
    
    ch4_indices = [c_idx] + ch3_h_indices + [diss_h_idx]
    
    return initial, ch4_indices


def prepare_neb_images(final_structure, n_images=11):
    """
    Prepare the NEB images with an initial structure generated from the final structure
    """
    initial_structure, ch4_indices = generate_initial_ch4_structure(final_structure)
    
    images = [initial_structure.copy()]
    for i in range(n_images - 2):
        images.append(initial_structure.copy())
    images.append(final_structure.copy())
    
    neb = NEB(images)
    neb.interpolate()
    
    return images, ch4_indices


def create_trajectory_visualization(images, ch4_indices, output_path):
    """
    Create a visualization of the reaction path by showing all CH4/CH3 positions
    overlaid on the final structure
    """
    final_structure = images[-1].copy()
    
    ch_positions = []
    ch_symbols = []
    
    for img in images:
        for idx in ch4_indices:
            ch_positions.append(img.positions[idx])
            ch_symbols.append(img.get_chemical_symbols()[idx])
    
    from ase import Atoms
    trajectory_vis = final_structure.copy()
    
    for pos, sym in zip(ch_positions, ch_symbols):
        trajectory_vis += Atoms(sym, positions=[pos])
    
    write_vasp(output_path, trajectory_vis, direct=True)
    
    return trajectory_vis


def main():
    base_input_folder = './neb_input_relaxed'
    output_folder = './neb_interpolated'
    os.makedirs(output_folder, exist_ok=True)
    
    o_numbers = [0]
    for o in o_numbers:
        indir = os.path.join(base_input_folder, f"NEB_Pd45O{o}")
        outdir = os.path.join(output_folder, f"NEB_Pd45O{o}")
        os.makedirs(outdir, exist_ok=True)
        
        site_nums = get_site_numbers(indir)
        
        for site_num in site_nums:
            fin_file = os.path.join(indir, f"Pd45O{o}_CH4_site{site_num}_final.vasp")
            traj_path = os.path.join(outdir, f"site{site_num}_neb_initial.traj")
            log_path = os.path.join(outdir, f"site{site_num}_neb.log")

            print(f"Processing: O{o}, site {site_num}")
            
            final_structure = read_vasp(fin_file)
            
            images, ch4_indices = prepare_neb_images(final_structure, n_images=11)
            
            initial_structure = images[0]
            write_vasp(os.path.join(outdir, f"site{site_num}_generated_initial.vasp"), initial_structure, direct=True)
            
            for i, img in enumerate(images):
                write_vasp(os.path.join(outdir, f"site{site_num}_image_{i:02d}.vasp"), img, direct=True)
            
            vis_path = os.path.join(outdir, f"site{site_num}_trajectory_visualization.vasp")
            create_trajectory_visualization(images, ch4_indices, vis_path)
            
            print(f"Created trajectory visualization: {vis_path}")


if __name__ == "__main__":
    main()
