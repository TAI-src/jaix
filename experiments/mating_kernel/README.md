# Mating Kernel

## Entrypoint

All experiments can be started via `run_exp.py`.

## Slurm experiments

### Initial setup

Create `pyproject.toml` in the folder.

```{toml}
[project]
name = "mating_kernel"
version = "0.1.0"
description = "Add your description here"
readme = "README.md"
requires-python = ">=3.12.13"
dependencies = [
    "cobi",
    "tai-jaix",
    "pydantic",
]

[tool.uv.sources]
cobi = { path = "../../deps/cobi" }
```

### Install / Update packages

```{bash}
module load gcc uv
uv lock --upgrade --python 3.12
uv sync --python 3.12
```

### Basic Slurm script

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
PYTHONUNBUFFERED=1 .venv/bin/python <fill_in_commands>
```

- The array id is available as `$SLURM_ARRAY_TASK_ID` in the script.
- The `%5`in the array specification imposes a limit of 5 concurrent jobs.

## Experiment types

```{bash}
python run_exp.py <exp_typ> ...
[--cobi] # include cobi problems
[--re] # include RE problems
[--num_objectives NUM_OBJECTIVES [NUM_OBJECTIVES ...]] # only include problems with the specified number of objectives
[--constrained] # include constrained problems
[--seed SEED] # starting seed for random number generator
[--reps REPS] # number of repetitions for each problem (and setting)
[--num_batches NUM_BATCHES] # number of batches to split the experiments into
--out_dir OUT_DIR # output directory for results
[--nth NTH [NTH ...]] # which of the num_batches to run (0-indexed)
```

### perf: Performance experiments

```{bash}
python run_exp.py perf ...
[--alg_name {NSGA2} [{NSGA2} ...]] # which algorithms to run
[--selector SELECTOR [SELECTOR ...]] # which selection operators to use
[--n_gen N_GEN] # number of generations
```

```{sh}
#!/bin/bash
#SBATCH -p medium96s
#SBATCH --job-name=perf_nsga2
#SBATCH --array=0-359%12
#SBATCH --cpus-per-task=4
#SBATCH -t 04:00:00
#SBATCH --output=logs/perf_nsga2/%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=thehedgeify@gmail.com

module load gcc uv
PYTHONUNBUFFERED=1 .venv/bin/python run_exp.py perf --nth "$SLURM_ARRAY_TASK_ID" --out_dir perf_nsga2 --cobi --re --num_objectives 2 --seed 1337 --reps 30 --num_batches 360 --n_gen 1000
```

We are running on all problems with 2 objectives, i.e. 7 cobi problems plus 5 RE problems, with 30 repetitions each, resulting in 360 tasks. Each task runs for 1000 generations.
