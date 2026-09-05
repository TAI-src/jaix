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

```{bash}
sbatch jobscript.sh
```
