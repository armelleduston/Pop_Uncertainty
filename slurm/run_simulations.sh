#!/bin/bash
#SBATCH --job-name=pop_sim
#SBATCH --output=logs/pop_sim_%A_%a.out
#SBATCH --error=logs/pop_sim_%A_%a.err
#SBATCH --array=0-63
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

# Load R if required by your cluster (optional)
module load R || true

cd "$SLURM_SUBMIT_DIR"

# Parameter grids (edit as needed)
etas=(0.1 0.5 1 2)
taus=(0.1 0.3 0.5 0.7)
kappas=(0.1 1 10 25)

E_LEN=${#etas[@]}
T_LEN=${#taus[@]}
K_LEN=${#kappas[@]}
TOTAL=$((E_LEN * T_LEN * K_LEN))

IDX=${SLURM_ARRAY_TASK_ID:-0}
if (( IDX >= TOTAL )); then
  echo "Task index $IDX >= TOTAL $TOTAL; adjust --array or parameter grids." >&2
  exit 1
fi

eta_index=$(( IDX / (T_LEN * K_LEN) ))
rem=$(( IDX % (T_LEN * K_LEN) ))
tau_index=$(( rem / K_LEN ))
kappa_index=$(( rem % K_LEN ))

eta=${etas[$eta_index]}
tau=${taus[$tau_index]}
kappa=${kappas[$kappa_index]}

mkdir -p logs results

echo "Starting simulation for tau=$tau kappa=$kappa eta=$eta (task $SLURM_ARRAY_TASK_ID)"
Rscript Simulation_file.R "$tau" "$kappa" "$eta"
echo "Finished simulation for tau=$tau kappa=$kappa eta=$eta"
