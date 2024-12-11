from jaix.utils import launch_jaix_experiment
import os
from wandb.sdk import launch
from typing import Dict, Optional
import shutil

xconfig = {
    "jaix.ExperimentConfig": {
        "env_config": {
            "jaix.EnvironmentConfig": {
                "suite_class": "jaix.suite.coco.COCOSuite",
                "suite_config": {
                    "jaix.suite.coco.COCOSuiteConfig": {
                        "env_config": {
                            "jaix.env.singular.ECEnvironmentConfig": {
                                "budget_multiplier": 1000,
                            },
                        },
                        "suite_name": "bbob",
                        "suite_instance": "instances: 1",
                        "suite_options": "function_indices: 1,2 dimensions: 2,3",
                        "num_batches": 1,
                        "current_batch": 0,
                        "output_folder": "test_run",
                    },
                },
                "env_wrappers": None,
                "comp_config": None,
                "seed": None,
            },
        },
        "runner_class": "jaix.runner.ask_tell.ATRunner",
        "runner_config": {
            "jaix.runner.ask_tell.ATRunnerConfig": {
                "max_evals": 4000,
                "disp_interval": 50,
            },
        },
        "opt_class": "jaix.runner.ask_tell.ATOptimiser",
        "opt_config": {
            "jaix.runner.ask_tell.ATOptimiserConfig": {
                "strategy_class": "jaix.runner.ask_tell.strategy.CMA",
                "strategy_config": {
                    "jaix.runner.ask_tell.strategy.CMAConfig": {
                        "sigma0": 2,
                    },
                },
                "init_pop_size": 1,
                "stop_after": 10000,
            },
        },
        "log_level": 20,
    },
}


def test_launch_jaix_experiment(config: Optional[Dict] = None):
    if not config:
        prev_mode = os.environ.get("WANDB_MODE", "online")
        os.environ["WANDB_MODE"] = "offline"
        run = launch_jaix_experiment(run_config=xconfig, project="ci-cd")
    else:
        run = launch_jaix_experiment(run_config=config)

    # Remove logging files
    shutil.rmtree(run.dir)

    if not config:
        # Iff available, reset the mode
        os.environ["WANDB_MODE"] = prev_mode


if __name__ == "__main__":
    # This is to test launch from wandb
    if not os.environ.get("WANDB_CONFIG", None):
        raise RuntimeError("Needs to be launched from wandb")
    run_config = launch.load_wandb_config()
    test_launch_jaix_experiment(run_config)
