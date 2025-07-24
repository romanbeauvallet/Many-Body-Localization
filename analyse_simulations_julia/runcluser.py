#!/usr/bin/env python3
import json
import subprocess
from pathlib import Path
import numpy as np

# =========================================== donner les valeurs de désorder

def sampling_disorder(x, h, y, sigma, H):
    """
    x -- point
    h -- upper boundary of the disorder magnitude
    y -- peak center index
    sigma -- peak width
    H -- amplitude of the peak

    Returns the value of the disorder distribution at x
    """
    if x < 0 or x > h:
        return 0.0
    return 1 + H * np.exp(-((x - y) ** 2) / (2 * sigma ** 2))


def rejection_sample(N, X, y, init, sigma=None, A=2.0):
    """
    N -- int: number of samples
    X -- float: maximum disorder magnitude
    y -- float: peak center
    init -- int: seed

    Optional:
    sigma -- float: width of the peak (default: 0.005 * X)
    A -- float: amplitude of the peak (default: 5.0)

    Returns:
    np.ndarray of shape (N,) with disorder samples
    """
    if sigma is None:
        sigma = 0.005 * X

    rng = np.random.default_rng(init)
    samples = []

    max_density = 1 + A  # upper bound of sampling_disorder

    while len(samples) < N:
        x = rng.uniform(0, X)
        u = rng.uniform(0, max_density)
        if u < sampling_disorder(x, X, y, sigma, A):
            samples.append(round(x, 4))

    return np.array(samples)

h = rejection_sample(20, 6, 2, 314159265, sigma=None, A=2.0)

# =========================================== faire le slurm et run

def generate_and_submit_jobs(seeds, jobname="sim_job", output_dir="slurm_jobs"):
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    Path("input").mkdir(exist_ok=True)

    for seed in seeds:
        slurm_script_path = output_path / f"{jobname}_seed{seed}.slurm"
        json_filename = f"input/input_seed{seed}.json"
        log_filename = f"{jobname}_seed{seed}.log"
        out_filename = f"{jobname}_seed{seed}.out"

        # Contenu du fichier JSON
        input_data = {
            "N": 100,
            "Ndisorder": 100,
            "J": 1,
            "D0": 10,
            "maximum value for beta": 10,
            "step value for beta list": 1,
            "site measure": 50,
            "cutoff": 1e-15,
            "gamma": 0,
            "max bond dimension": 300,
            "Trotter-Suzuki step": 1e-3,
            "fixed disorder": 0,
            "maximum disorder magnitude": 5,
            "nsweep range": [100, 1000, 100],
            "phase transition point": 2,
            "fixed seed": seed,
            "fixed number of sweep": 3000,
            "initialization": "neel",
            "gammescale": 0.8,
            "axis": "z",
            "gammelength": [50, 100, 10],
            "liste step Trotter-Suzuki": [1e-1, 1e-2, 1e-3, 1e-4, 1e-5],
            "number sweep dmrg": 20,
            "noise": 1e-8,
            "model choisi": "XY",
            "beta fixe": 40,
            "savefile": "data/ErrorAgainstTau.json"
        }

        # Sauvegarde du fichier JSON
        with open(json_filename, "w") as f_json:
            json.dump(input_data, f_json, indent=4)

        # Génération du script SLURM
        slurm_text = f"""#!/bin/bash
#SBATCH --job-name={jobname}_{seed}
#SBATCH --output={out_filename}
#SBATCH --error={log_filename}
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu
#SBATCH --ntasks=1

INPUT={json_filename}

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=false

echo STARTING AT $(date '+%A %d %B %Y %H:%M')
echo HOST $(hostname)
echo JOB ID $SLURM_JOB_ID

if ! jq empty $INPUT; then
  echo "ERROR: $INPUT is not valid"
  exit 1
fi

module load modules/2.3-20240529 git openblas/threaded-0.3.26

julia -O2 --project=/mnt/home/rbeauvallet/Many-Body-Localization --startup-file=no --math-mode=fast /mnt/home/rbeauvallet/Many-Body-Localization/run/Scaling_for_cluster.jl $INPUT

echo FINISHED AT $(date '+%A %d %B %Y %H:%M')
"""

        with open(slurm_script_path, "w") as f_slurm:
            f_slurm.write(slurm_text)

        # Envoi du job
        print(f"Submitting job for seed = {seed}")
        subprocess.run(["sbatch", str(slurm_script_path)])


# Exemple d’utilisation
if __name__ == "__main__":
    seed_list = [12345, 23456, 34567, 45678]  # Modifie cette liste selon tes besoins
    generate_and_submit_jobs(seed_list)
