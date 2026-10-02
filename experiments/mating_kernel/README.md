# Mating Kernel

## Entrypoint

All experiments can be started via `run_exp.py`.

## Slurm experiments

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
