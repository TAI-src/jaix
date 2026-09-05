# Moomap experiments

## Slurm script

```{sh}
#!/bin/bash
#SBATCH -p standard96s
#SBATCH --job-name=nsga3
#SBATCH --array=0-22%5
#SBATCH --cpus-per-task=4
#SBATCH -t 03:00:00
#SBATCH --output=logs/nsga3_%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>

module load gcc uv
uv run nsga3_experiment.py --num_independent_runs 1 --num_generations 1000 --problem_idx "$SLURM_ARRAY_TASK_ID"
```

To start the jobs, run the following. Done in 2 batches to set different time limits for the second batch of jobs since they take longer.

```{bash}
for i in {1..30}; do
    sbatch --array=0-11,13-20%5 jobscript.sh
    sbatch --array=12,21-22%5 --time=05:00:00 jobscript.sh
done
```

For retries

```{bash}
declare -A retries=(
    [12]=8
    [21]=12
    [22]=30
)

for task_id in "${!retries[@]}"; do
    for ((i=0; i<retries[$task_id]; i++)); do
        sbatch --array="$task_id" --time=05:00:00 jobscript.sh
    done
done
```

## Update package and repo on cluster

```{bash}
cd $PROJECT
cd jaix/experiments/moomap
git pull

module load gcc uv
uv lock --upgrade
uv sync
```
