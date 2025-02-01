# Ensemble selection

Corresponding singular environment [HPOEnvironment](../../jaix/env/singular/hpo_env.py) and composite environment [SwitchingEnvironment](../../jaix/env/composite/switching_environment.py) where relevant.

## Configuration Details

* [singular config](./binary.json)
* [warm start config](./binary_warm.json)
* [composite config](./binary_comp.json) (not in paper)

Important details:
* `task_type=C1`: Corresponds to binary classification
* `warm_start_strategy`: None or best, depending on experiment
