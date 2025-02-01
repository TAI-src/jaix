# Brachytherapy / RBF experiments

Corresponding singular environment [ECEnvironment](../../jaix/env/singular/ec_env.py) using problem [RBFFit](../../jaix/env/utils/problem/rbf_fit.py).

## Configuration Details
Full [config](brachy.json). Important details highlighted below:
* search space dimension between 15-25, depending on the instance
* `num_measure_points`: number of measurement points for value passed to optimiser
* `num_true_measurement_points`: number of measurement points for logging
* `noisy = true`, which means that the measurement points are different(resampled) for the values passed to the optimiser vs the logging
* `precision= 1e-8`: Environment terminates when the error between target function and RBF is less than `1e-8`
* `max_evals=5000`: Maximum number of evaluations

## Details

see paper

