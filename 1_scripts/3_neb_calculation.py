"""
Perform NEB (Nudged Elastic Band) calculations for CH₄ activation on Pd–CeO₂ slabs
with varying oxygen content.

The script:
1. Loads initial and final states based on pre-generated VASP images.
2. Optimizes the initial and final geometries to ensure physical feasibility.
3. Runs NEB calculations using the pre-optimized images, with force convergence criteria.
4. Saves the optimized NEB images and energy profiles for each site.
5. Generates and saves plots of energy vs. image number and max force vs. image number.
6. Identifies and saves the transition state (highest energy image) for each run.

Dependencies: ASE, NumPy, SciPy.

Author: Seokhyun Choung
Licence: MIT
"""


import os
import re
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import read, write
from ase.optimize import LBFGS, FIRE
from ase.neb import NEB
from ase.io.vasp import read_vasp
from ase.calculators.singlepoint import SinglePointCalculator
import glob
import shutil
import sys
import signal
import traceback

# Base output directories
base_output_dir = './neb_output'
best5_output_dir = './neb_output_all'
os.makedirs(base_output_dir, exist_ok=True)
os.makedirs(best5_output_dir, exist_ok=True)

class GracefulInterruptHandler:
    def __init__(self):
        self.interrupted = False
        self.released = False
        self.original_handler = None

    def __enter__(self):
        self.original_handler = signal.getsignal(signal.SIGINT)
        signal.signal(signal.SIGINT, self.handler)
        return self

    def handler(self, sig, frame):
        self.interrupted = True
        print("\nInterrupt received, finishing current operation gracefully...")
        signal.signal(signal.SIGINT, self.original_handler)

    def __exit__(self, type, value, traceback):
        if self.interrupted:
            print("Exiting due to keyboard interrupt")
        signal.signal(signal.SIGINT, self.original_handler)
        self.released = True

def log_message(message, logfile):
    """Log message to file and print"""
    with open(logfile, "a") as log_file:
        log_file.write(message + "\n")
    print(message)

def get_opt_energy(atoms, fmax=0.05, steps=100, opt_mode="normal", traj_file=None, logfile=None):
    """Optimize structure"""
    log_message(f"Optimizing structure with fmax={fmax}, steps={steps}, mode={opt_mode}...", logfile)
    atoms.set_calculator(calculator)
    opt = LBFGS(atoms, trajectory=traj_file) if opt_mode == "normal" else LBFGS(atoms)
    
    with GracefulInterruptHandler() as h:
        try:
            opt.run(fmax=fmax, steps=steps)
            log_message("Optimization complete.", logfile)
        except Exception as e:
            if h.interrupted:
                log_message("Optimization interrupted by user.", logfile)
            else:
                log_message(f"Error during optimization: {str(e)}", logfile)
                traceback.print_exc()
    
    return atoms

def find_best_sites(o_folder, o_number, n_best=100):
    """Find best sites based on final structure energy"""
    site_numbers = get_site_numbers(o_folder)
    temp_dir = os.path.join(base_output_dir, f"temp_energy_calc_O{o_number}")
    os.makedirs(temp_dir, exist_ok=True)
    temp_log = os.path.join(temp_dir, "energy_calc.log")
    
    site_energies = {}
    for site_number in site_numbers:
        try:
            final_file = os.path.join(o_folder, f"site{site_number}_image_10.vasp")
            if not os.path.exists(final_file):
                continue
            final_atoms = read_vasp(final_file)
            final_atoms.calc = calculator
            site_energies[site_number] = final_atoms.get_potential_energy()
        except Exception as e:
            print(f"Error calculating energy for site {site_number}: {str(e)}")
    
    sorted_sites = sorted(site_energies.items(), key=lambda x: x[1])
    best_sites = [site for site, _ in sorted_sites[:n_best]]
    return best_sites

def run_neb_with_vasp_images(images_dir, site_number, o_number, n_images=11, relax_endpoints=True, fmax_neb=0.05, steps_neb=500, logfile=None):
    """Run NEB with pre-generated VASP images"""
    output_dir = os.path.dirname(logfile)
    images = []
    
    for i in range(n_images):
        image_file = os.path.join(images_dir, f"site{site_number}_image_{i:02d}.vasp")
        if os.path.exists(image_file):
            img = read_vasp(image_file)
            images.append(img)
    
    if relax_endpoints:
        images[0] = get_opt_energy(images[0], fmax=0.01, steps=200, logfile=logfile)
        images[-1] = get_opt_energy(images[-1], fmax=0.01, steps=200, logfile=logfile)
    
    for img in images:
        img.calc = calculator
    
    neb = NEB(images, k=0.03, climb=True, allow_shared_calculator=True)
    optimizer = FIRE(neb, logfile=logfile)
    with GracefulInterruptHandler() as h:
        try:
            optimizer.run(fmax=fmax_neb, steps=steps_neb)
            log_message("NEB optimization complete.", logfile)
        except Exception as e:
            if h.interrupted:
                log_message("NEB optimization interrupted by user.", logfile)
            else:
                log_message(f"Error during NEB optimization: {str(e)}", logfile)
                traceback.print_exc()
    
    return images

def plot_neb(images, output_dir, o_number, site_number, logfile=None):
    """Plot energy and max force for each image"""
    forces = [calc_max_force(img) for img in images]
    energies = [img.get_potential_energy() for img in images]
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.set_xlabel('Image')
    ax1.set_ylabel('Energy (eV)', color='tab:blue')
    ax1.plot(range(len(energies)), energies, marker='o', color='tab:blue')
    ax2 = ax1.twinx()
    ax2.set_ylabel('Max Force (eV/Å)', color='tab:red')
    ax2.plot(range(len(forces)), forces, marker='x', linestyle='--', color='tab:red')
    plt.title(f"NEB Results: O{o_number} Site {site_number}")
    plot_path = os.path.join(output_dir, f"O{o_number}_site{site_number}_neb_plot.png")
    plt.savefig(plot_path)
    plt.close()
    return energies, forces

def save_transition_state(images, output_dir, o_number, site_number, logfile=None):
    """Save transition state (highest energy image)"""
    energies = [img.get_potential_energy() for img in images]
    ts_index = np.argmax(energies)
    ts_image = images[ts_index]
    transition_state_file = os.path.join(output_dir, f"O{o_number}_site{site_number}_ts.vasp")
    write(transition_state_file, ts_image)
    energy_file = os.path.join(output_dir, f"O{o_number}_site{site_number}_energies.txt")
    with open(energy_file, 'w') as f:
        f.write("# NEB Energy Profile\n")
        for i, e in enumerate(energies):
            f.write(f"{i:6d}  {e:12.6f}  {e-energies[0]:12.6f}\n")
    return ts_index, energies[ts_index]

def create_output_directory(o_number, site_number):
    """Create output directory for each O and site combination"""
    site_dir = os.path.join(best5_output_dir, f"PdO{o_number}_site{site_number}")
    os.makedirs(site_dir, exist_ok=True)
    return site_dir

def get_site_numbers(folder):
    """Get site numbers from folder based on image files"""
    files = os.listdir(folder)
    site_numbers = set()
    for fname in files:
        match = re.search(r'site(\d+)_image_\d+\.vasp', fname)
        if match:
            site_numbers.add(int(match.group(1)))
    return sorted(site_numbers)

def create_combined_plot(o_number, site_numbers, output_dir):
    """Create combined plot of energy profiles for all sites"""
    plt.figure(figsize=(12, 8))
    for site in site_numbers:
        try:
            energy_file = os.path.join(output_dir, f"PdO{o_number}_site{site}", f"O{o_number}_site{site}_energies.txt")
            if os.path.exists(energy_file):
                data = np.loadtxt(energy_file, skiprows=2)
                images = data[:, 0]
                rel_energies = data[:, 2]
                plt.plot(images, rel_energies, 'o-', label=f'Site {site}')
        except Exception as e:
            print(f"Error loading data for site {site}: {str(e)}")
    
    plt.xlabel('Image')
    plt.ylabel('Relative Energy (eV)')
    plt.title(f'NEB Energy Profiles for O{o_number}')
    plt.grid(True)
    plt.legend()
    plt.savefig(os.path.join(output_dir, f"O{o_number}_combined_profiles.png"))
    plt.close()

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nProcess interrupted by user. Exiting gracefully.")
        sys.exit(1)
