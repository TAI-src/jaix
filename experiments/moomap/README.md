# Moomap experiments

## Slurm experiments

Basic slurm script

```{sh}
#!/bin/bash
#SBATCH -p standard96s
#SBATCH --job-name=nsga3
#SBATCH --array=0-22%5
#SBATCH --cpus-per-task=4
#SBATCH -t 02:00:00
#SBATCH --output=logs/nsga3_%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>

module load gcc uv
PYTHONUNBUFFERED=1 uv run <fill_in_commands>
```

Instead of `uv run`, it is also possible to directly use the venv that is created in advance (e.g. in the login node). `.venv/bin/python <fill_in_commands>`. This might be more stable on a compute node without internet access but requires that the venv is created in advance and that the correct python version is used.

To start 30 batch jobs (with 23 tasks each), run the following. The SBATCH values can be overridden as shown.

- The array id is available as `$SLURM_ARRAY_TASK_ID` in the script.
- The `%5`in the array specification imposes a limit of 5 concurrent jobs.

```{bash}
for i in {1..30}; do
  sbatch --array=0-11,13-20%5 jobscript.sh
  sbatch --array=12,21-22%5 --time=05:00:00 jobscript.sh
done
```

In case specific task ids need to be rerun:

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

Show 20 most recent completed jobs with their state:

```{bash}
sacct -u $USER -X --starttime 1970-01-01 --format=JobID,JobName,State,Elapsed,ExitCode | tail -20
```

### Offspring generation experiments

```{bash}
PYTHONUNBUFFERED=1 uv run nsga3_experiment.py --num_independent_runs 1 --num_generations 1000 --out_dir offspring_results --problem_idx "$SLURM_ARRAY_TASK_ID"
```

There are 23 possible task ids (0-22) for the problem index. Problem 0-6 are the cobi problems, 7-22 are the REProblems (unconstrained). The following table shows the recommended execution times for each task id.

- 22: more than 17:00:00
- 21: 10:00:00 mostly fine
- 12: 07:00:00 is enough
- rest: 02:00:00 is enough

### Postprocess data

```{bash}
PYTHONUNBUFFERED=1 uv run postprocess.py --results_dir offspring_results --out_dir ppdata --skip_plots --problem_idx "$SLURM_ARRAY_TASK_ID"
```

There are 23 possible task ids (0-22) for the problem index. Problem 0-6 are the cobi problems, 7-22 are the REProblems (unconstrained).

10 minutes are enough.

### Prediction feature importance experiments

```{bash}
PYTHONUNBUFFERED=1 uv run pred.py --batch_id $SLURM_ARRAY_TASK_ID --data_dir ppdata --output_dir pred_res
```

By default, there are 32 scenarios (input and target column combinations) and one data file per problem. This can be restricted by specifying `--scenario_ids` and `--file_ids`. The number of available batches is the product of the number of scenarios and the number of files. So for the full spread of 32 scenarios and 23 files, there are 736 batches.

Recommended batch time: 30:00 should be fine

## Update package and repo on cluster

```{bash}
cd $PROJECT
cd jaix/experiments/moomap
git pull

module load gcc uv
uv lock --upgrade --python 3.12
uv sync --python 3.12
```
