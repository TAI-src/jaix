from jaix.utils import launch_jaix_experiment, wandb_logger, wandb_init
import os
import shutil
from ttex.config import ConfigFactory as CF
from copy import deepcopy

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
        "logging_config": {
            "jaix.LoggingConfig": {
                "log_level": 10,
            }
        },
    },
}


def test_wandb_logger():
    exp_config = CF.from_dict(xconfig)
    nexp_config = wandb_logger(exp_config, "dummy_run", "dummy_name")
    assert "dummy_name" in nexp_config.logging_config.dict_config["loggers"]
    logger_config = nexp_config.logging_config.dict_config["loggers"]["dummy_name"]
    assert "wandb_handler" in logger_config["handlers"]
    assert "wandb_handler" in nexp_config.logging_config.dict_config["handlers"]
    logging_wrapper_tuple = nexp_config.env_config.env_wrappers[-1]
    assert logging_wrapper_tuple[1].logger_name == "dummy_name"


def test_wandb_init():
    prev_mode = os.environ.get("WANDB_MODE", "online")
    os.environ["WANDB_MODE"] = "offline"
    run = wandb_init(run_config=deepcopy(xconfig), project="ci-cd")
    assert run.mode == "dryrun"
    shutil.rmtree(run.dir)
    run.finish()

    os.environ["WANDB_MODE"] = prev_mode


def test_launch_jaix_experiment():
    prev_mode = os.environ.get("WANDB_MODE", "online")
    os.environ["WANDB_MODE"] = "offline"
    data_dir, exit_code = launch_jaix_experiment(
        run_config=deepcopy(xconfig), project="ci-cd"
    )

    # Remove logging files
    shutil.rmtree(data_dir)
    os.environ["WANDB_MODE"] = prev_mode

    assert exit_code == 0
