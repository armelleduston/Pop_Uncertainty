#!/bin/bash
#SBATCH --job-name=combined_sim
#SBATCH --output=slurm/logs/combined_sim_%A_%a.out
#SBATCH --error=slurm/logs/combined_sim_%A_%a.err
#SBATCH --array=1-10                # 10 tasks (adjust as needed)
#SBATCH --ntasks=1                  # 1 task per node
#SBATCH --cpus-per-task=20          # cores per node (adjust based on cluster)
#SBATCH --mem=32G                   # memory per node (adjust as needed)
#SBATCH --time=02:00:00             # max 2 hours per task
#SBATCH --partition=standard        # adjust to your cluster's partition name

# Load R module (adjust module name based on your cluster)
module load R/4.3.0

# Set R library path to use custom package installation
export R_LIBS_USER=~/R_libs/pop_uncertainty

# Print job information
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running on node: $(hostname)"
echo "Number of cores: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Create logs directory if it doesn't exist
mkdir -p slurm/logs

# Run the R script
Rscript full_sim/combined_simulation.R

echo "End time: $(date)"
echo "Task $SLURM_ARRAY_TASK_ID completed"
