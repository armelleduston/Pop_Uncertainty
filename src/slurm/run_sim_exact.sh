#!/bin/bash
#SBATCH --job-name=simulations
#SBATCH --output=sbatch_logs/run_sim_%A_%a.out
#SBATCH --error=sbatch_logs/run_sim_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=1-100%10

# Project directory (defaults to the directory where sbatch is submitted)
PROJECT_DIR=${SLURM_SUBMIT_DIR:-$PWD}

# Ensure folders exist
mkdir -p "$PROJECT_DIR/sbatch_logs"
mkdir -p "$PROJECT_DIR/results"

# Load R module if your cluster uses modules (edit as needed)
module load R || true

cd "$PROJECT_DIR"

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
echo "Starting simulation task ${TASK_ID} on $(hostname) at $(date)"

# Run the R script with a reproducible seed and save the samples object
srun Rscript -e "set.seed(${TASK_ID}); source('Simulation_file.R'); dir.create('results', showWarnings=FALSE); saveRDS(samples2, file=paste0('results/samples_', ${TASK_ID}, '.rds'))"

echo "Finished simulation task ${TASK_ID} at $(date)"