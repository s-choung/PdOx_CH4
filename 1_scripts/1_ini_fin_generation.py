"""Generate initial and final geometries for CH₄ activation NEB runs on
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

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, List, Sequence

import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.io import read, write
from ase.neighborlist import neighbor_list
from scipy.spatial.transform import Rotation as R

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ──────────────────────────────────────────────────────────────────────────────
BEST_STRUCTURE_DIR = Path(
    "input_slab_path"
)
OUTPUT_DIR = Path(
    "neb_input"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

O_NUMS: Sequence[int] = (0, 20, 45)          # slab oxidation states
CLEARANCE_CH4: float = 3.0                   # Å C–Pd distance in initial state
CLEARANCE_CH3: float = 2.0                   # Å C–Pd distance in final state
SURFACE_THRESHOLD_Z: float = 6.0             # Å Pd considered “surface Pd”

# Distance cut‐offs for overlap test (Å)
H_MIN_DIST: float = 0.8
DEFAULT_MIN_DIST: float = 1.8

logging.basicConfig(
    format="[%(levelname)s] %(message)s", level=logging.INFO, force=True
)

# ──────────────────────────────────────────────────────────────────────────────
# UTILITY FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

def load_slab(o_num: int) -> Atoms | None:
    """Return ASE Atoms of the slab with *o_num* oxygens, or *None* if missing."""
    poscar = BEST_STRUCTURE_DIR / f"Pd45_O{o_num}" / "POSCAR.vasp"
    if not poscar.is_file():
        logging.warning("POSCAR not found for O%s", o_num)
        return None
    return read(poscar)


def surface_pd_indices(atoms: Atoms) -> List[int]:
    """Indices of Pd atoms within *SURFACE_THRESHOLD_Z* of the highest Pd."""
    zs = atoms.positions[:, 2]
    pd_idx = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "Pd"]
    if not pd_idx:
        return []
    max_z = zs[pd_idx].max()
    return [i for i in pd_idx if (max_z - zs[i]) <= SURFACE_THRESHOLD_Z]


def pd_centroid(atoms: Atoms, pd_idx: Sequence[int]) -> np.ndarray:
    """Centroid used as upward reference: (mean x,y, min z) of Pd atoms."""
    pd_pos = atoms.positions[pd_idx]
    centroid_xy = pd_pos[:, :2].mean(axis=0)
    min_z = pd_pos[:, 2].min()
    return np.array([*centroid_xy, min_z])


def _align_fragment(fragment: Atoms, up_vector: np.ndarray) -> None:
    """Rotate *fragment* so that the mean (fragment-H) vector aligns with *up_vector*."""
    c_pos = fragment.positions[0]
    current = (fragment.positions[1:] - c_pos).mean(axis=0)
    current /= np.linalg.norm(current)
    up = up_vector / np.linalg.norm(up_vector)
    axis = np.cross(current, up)
    angle = np.arccos(np.clip(np.dot(current, up), -1.0, 1.0))
    if np.linalg.norm(axis) > 1e-3 and angle > 1e-3:
        rot = R.from_rotvec(axis / np.linalg.norm(axis) * angle)
        fragment.positions[:] = rot.apply(fragment.positions - c_pos) + c_pos


def make_ch3(up_vector: np.ndarray) -> Atoms:
    """Return a CH₃ fragment pointing opposite *up_vector* (C index 0)."""
    ch4 = molecule("CH4")
    c_pos = ch4.positions[0]
    h_vectors = ch4.positions[1:] - c_pos
    # Remove the H that points *most* towards −up_vector to make CH3
    align = [np.dot(h / np.linalg.norm(h), -up_vector / np.linalg.norm(up_vector)) for h in h_vectors]
    remove_idx = 1 + int(np.argmax(align))
    ch3 = ch4[[i for i in range(len(ch4)) if i != remove_idx]]
    _align_fragment(ch3, up_vector)
    return ch3


def make_ch4(up_vector: np.ndarray) -> Atoms:
    ch4 = molecule("CH4")
    _align_fragment(ch4, up_vector)
    return ch4


def place_fragment(slab: Atoms, frag: Atoms, pd_pos: np.ndarray, up_vector: np.ndarray, clearance: float) -> Atoms:
    """Return *slab + frag* with C sitting *clearance* Å above *pd_pos*."""
    shift = pd_pos + up_vector * clearance - frag.positions[0]
    frag.positions += shift
    return slab + frag


def min_allowed_dist(sym1: str, sym2: str) -> float:
    return H_MIN_DIST if "H" in (sym1, sym2) else DEFAULT_MIN_DIST


def has_bad_overlap(ref: Atoms, probe: Atoms) -> bool:
    ref_pos, probe_pos = ref.positions, probe.positions
    ref_sym, probe_sym = ref.get_chemical_symbols(), probe.get_chemical_symbols()
    for i, p in enumerate(probe_pos):
        for j, r in enumerate(ref_pos):
            if np.linalg.norm(p - r) < min_allowed_dist(probe_sym[i], ref_sym[j]):
                return True
    return False


def nearest_O_or_Pd(atoms: Atoms, idx: int, exclude: Iterable[int] = (), cutoff: float = 3.5) -> int | None:
    sym = atoms.get_chemical_symbols()
    i_list, j_list, dist = neighbor_list("ijd", atoms, cutoff=cutoff)
    neigh = [(j, d) for i, j, d in zip(i_list, j_list, dist) if i == idx and j not in exclude]
    o_neigh = [(j, d) for j, d in neigh if sym[j] == "O"]
    pd_neigh = [(j, d) for j, d in neigh if sym[j] == "Pd"]
    if o_neigh:
        return min(o_neigh, key=lambda x: x[1])[0]
    if pd_neigh:
        return min(pd_neigh, key=lambda x: x[1])[0]
    return None


def place_H(atoms: Atoms, target_idx: int, centroid: np.ndarray, distance: float = 1.0) -> Atoms:
    target_pos = atoms.positions[target_idx]
    up_vector = (target_pos - centroid) / np.linalg.norm(target_pos - centroid)
    h_pos = target_pos + up_vector * distance
    return atoms + Atoms("H", positions=[h_pos])

# ──────────────────────────────────────────────────────────────────────────────
# MAIN ROUTINE
# ──────────────────────────────────────────────────────────────────────────────

def generate_neb_files() -> None:
    for o_num in O_NUMS:
        slab = load_slab(o_num)
        if slab is None:
            continue

        pd_surf = surface_pd_indices(slab)
        if not pd_surf:
            logging.warning("No surface Pd found for O%s", o_num)
            continue
        centroid = pd_centroid(slab, pd_surf)

        save_root = OUTPUT_DIR / f"NEB_Pd45O{o_num}"
        save_root.mkdir(exist_ok=True)

        logging.info("O%s: %d surface Pd atoms", o_num, len(pd_surf))

        for site_num, pd_idx in enumerate(pd_surf, start=1):
            pd_pos = slab.positions[pd_idx]
            up_vec = (pd_pos - centroid) / np.linalg.norm(pd_pos - centroid)

            # Build initial state (CH4)
            ch4 = make_ch4(up_vec)
            init_atoms = place_fragment(slab, ch4, pd_pos, up_vec, CLEARANCE_CH4)

            if has_bad_overlap(slab, init_atoms[len(slab):]):
                logging.debug("Site %d skipped (CH4 overlap)", pd_idx)
                continue

            # Build final state (CH3 + H)
            ch3 = make_ch3(up_vec)
            final_atoms = place_fragment(slab, ch3, pd_pos, up_vec, CLEARANCE_CH3)
            ch3_C_idx = len(slab)  # first atom of fragment in combined system
            target_idx = nearest_O_or_Pd(final_atoms, ch3_C_idx, exclude=[pd_idx])
            if target_idx is None:
                logging.debug("Site %d skipped (no target for H)", pd_idx)
                continue
            final_atoms = place_H(final_atoms, target_idx, centroid, distance=1.0)

            if has_bad_overlap(final_atoms[: len(slab)], final_atoms[len(slab) :]):
                logging.debug("Site %d skipped (final overlap)", pd_idx)
                continue

            # Write VASP files
            stem = f"Pd45O{o_num}_CH4_site{site_num}"
            write(save_root / f"{stem}_initial.vasp", init_atoms, vasp5=True, direct=True)
            write(save_root / f"{stem}_final.vasp", final_atoms, vasp5=True, direct=True)
            logging.info("Saved site %d (O%s): initial & final", site_num, o_num)


if __name__ == "__main__":
    generate_neb_files()
