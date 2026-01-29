# Combined Simulation for SLURM Cluster

## Overview
This setup runs 500 simulations across multiple nodes on a SLURM cluster. Each simulation runs both bench and no_bench models with 3 chains.

## Files Created
- `full_sim/combined_simulation.R` - Modified to work with SLURM array jobs
- `full_sim/setup_r_environment.R` - Installs all required R packages
- `slurm/submit_combined_simulation.sh` - SLURM submission script
- `full_sim/combine_results.R` - Script to combine results from all tasks

## How to Use

### 0. First-time Setup: Install R Packages

On the cluster, run once to install all required packages:
```bash
# Load R module
module load R/4.3.0

# Run setup script (takes 10-20 minutes)
Rscript full_sim/setup_r_environment.R
```

This installs packages to `~/R_libs/pop_uncertainty` which the SLURM script will use automatically.

### 1. Adjust SLURM Parameters
Edit `slurm/submit_combined_simulation.sh` and adjust these settings for your cluster:
- `#SBATCH --array=1-10` - Number of tasks (10 tasks = 50 runs each)
- `#SBATCH --cpus-per-task=20` - Cores per node
- `#SBATCH --mem=32G` - Memory per node
- `#SBATCH --partition=standard` - Your cluster's partition name
- `module load R/4.3.0` - Your cluster's R module

### 2. Create logs directory
```bash
mkdir -p slurm/logs
```

### 3. Submit the job
```bash
cd /path/to/Pop_Uncertainty
sbatch slurm/submit_combined_simulation.sh
```

### 4. Monitor progress
```bash
# Check job status
squeue -u $USER

# View output from a specific task
tail -f slurm/logs/combined_sim_JOBID_TASKID.out

# Check for errors
tail -f slurm/logs/combined_sim_JOBID_TASKID.err
```

### 5. After all tasks complete, combine results
```bash
Rscript full_sim/combine_results.R
```

This creates:
- `full_sim/results/no_bench_results_combined.csv`
- `full_sim/results/bench_results_combined.csv`

## Performance Notes
- With 10 tasks: Each runs ~50 simulations on one node
- With 20 cores per task: Uses `mclapply` to parallelize within each node
- Total time: Depends on cores/node, typically 30-90 minutes

## Adjusting Number of Tasks
To change how work is distributed:
- More tasks (e.g., `--array=1-20`): Faster overall, but more overhead
- Fewer tasks (e.g., `--array=1-5`): Longer runtime, less overhead
- Tasks automatically divide the 500 runs evenly
