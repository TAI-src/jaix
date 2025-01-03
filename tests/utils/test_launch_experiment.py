from jaix.utils import launch_jaix_experiment, wandb_logger, wandb_init
import os
import shutil
from ttex.config import ConfigFactory as CF
from copy import deepcopy
import pytest
import itertools


def get_config(suite="COCO", comp=False):
    xconfig = {
        "jaix.ExperimentConfig": {
            "runner_class": "jaix.runner.ask_tell.ATRunner",
            "runner_config": {
                "jaix.runner.ask_tell.ATRunnerConfig": {
                    "max_evals": 100,
                    "disp_interval": 50,
                },
            },
            "logging_config": {
                "jaix.LoggingConfig": {
                    "log_level": 10,
                }
            },
        },
    }
    if suite == "COCO":
        xconfig["jaix.ExperimentConfig"]["env_config"] = {
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
        }
    elif suite == "RBF":
        xconfig["jaix.ExperimentConfig"]["env_config"] = {
            "jaix.EnvironmentConfig": {
                "suite_class": "jaix.suite.ECSuite",
                "suite_config": {
                    "jaix.suite.ECSuiteConfig": {
                        "func_class": "jaix.env.utils.problem.RBFFit",
                        "func_config": {
                            "jaix.env.utils.problem.RBFFitConfig": {
                                "rbf_config": {
                                    "jaix.env.utils.problem.rbf.RBFAdapterConfig": {},
                                },
                                "precision": 1e-8,
                            },
                        },
                        "env_config": {
                            "jaix.env.singular.ECEnvironmentConfig": {
                                "budget_multiplier": 2,
                            },
                        },
                        "num_instances": 5,
                        "num_agg_instances": 2,
                    },
                },
                "env_wrappers": None,
                "comp_config": None,
                "seed": None,
            },
        }
    elif suite == "HPO":
        xconfig["jaix.ExperimentConfig"]["env_config"] = {
            "jaix.EnvironmentConfig": {
                "suite_class": "jaix.suite.Suite",
                "suite_config": {
                    "jaix.suite.SuiteConfig": {
                        "env_class": "jaix.env.singular.HPOEnvironment",
                        "env_config": {
                            "jaix.env.singular.HPOEnvironmentConfig": {
                                "training_budget": 100,
                                "task_type": "jaix.env.utils.hpo.TaskType.C1",
                                "repo_name": "D244_F3_C1530_3",
                                "load_predictions": False,
                                "cache": True,
                            },
                        },
                        "num_instances": 2,
                        "num_agg_instances": 3,
                    },
                },
                "env_wrappers": None,
                "comp_config": None,
                "seed": None,
            },
        }
    elif suite == "MMind":
        xconfig["jaix.ExperimentConfig"]["env_config"] = {
            "jaix.EnvironmentConfig": {
                "suite_class": "jaix.suite.Suite",
                "suite_config": {
                    "jaix.suite.SuiteConfig": {
                        "env_class": "jaix.env.singular.MastermindEnvironment",
                        "env_config": {
                            "jaix.env.singular.MastermindEnvironmentConfig": {
                                "num_slots_range": (4, 6),
                                "num_colours_range": (2, 3),
                                "sequential": False,
                                "max_guesses": 10,
                            },
                        },
                        "num_instances": 2,
                        "num_agg_instances": 3,
                    },
                },
                "env_wrappers": None,
                "comp_config": None,
                "seed": None,
            },
        }
    if comp:
        xconfig["jaix.ExperimentConfig"]["env_config"]["jaix.EnvironmentConfig"][
            "comp_config"
        ] = {
            "jaix.CompositeEnvironmentConfig": {
                "agg_type": "jaix.suite.AggType.INST",
                "comp_env_class": "jaix.env.composite.SwitchingEnvironment",
                "comp_env_config": {
                    "jaix.env.composite.SwitchingEnvironmentConfig": {
                        "switching_pattern_class": "jaix.env.utils.switching_pattern.SeqRegSwitchingPattern",
                        "switching_pattern_config": {
                            "jaix.env.utils.switching_pattern.SeqRegSwitchingPatternConfig": {
                                "wait_period": 100,
                            },
                        },
                        "real_time": False,
                    },
                },
            }
        }
        xconfig["jaix.ExperimentConfig"][
            "opt_class"
        ] = "jaix.runner.ask_tell.ATOptimiser"
        xconfig["jaix.ExperimentConfig"]["opt_config"] = {
            "jaix.runner.ask_tell.ATOptimiserConfig": {
                "strategy_class": "jaix.runner.ask_tell.strategy.ATBandit",
                "strategy_config": {
                    "jaix.runner.ask_tell.strategy.ATBanditConfig": {
                        "bandit_config": {
                            "jaix.runner.ask_tell.strategy.utils.BanditConfig": {
                                "epsilon": 0.1,
                                "min_tries": 10,
                                "exploit_strategy": "jaix.runner.ask_tell.strategy.utils.BanditExploitStrategy.MAX",
                            },
                        },
                        "opt_confs": [
                            {
                                "jaix.runner.ask_tell.ATOptimiserConfig": {
                                    "strategy_class": "jaix.runner.ask_tell.strategy.CMA",
                                    "strategy_config": {
                                        "jaix.runner.ask_tell.strategy.CMAConfig": {
                                            "sigma0": 2,
                                        },
                                    },
                                    "init_pop_size": 1,
                                    "stop_after": 400,
                                }
                            },
                            {
                                "jaix.runner.ask_tell.ATOptimiserConfig": {
                                    "strategy_class": "jaix.runner.ask_tell.strategy.CMA",
                                    "strategy_config": {
                                        "jaix.runner.ask_tell.strategy.CMAConfig": {
                                            "sigma0": 2,
                                        },
                                    },
                                    "init_pop_size": 1,
                                    "stop_after": 400,
                                }
                            },
                        ],
                    },
                },
                "init_pop_size": 1,
            },
        }
    else:
        xconfig["jaix.ExperimentConfig"][
            "opt_class"
        ] = "jaix.runner.ask_tell.ATOptimiser"
        xconfig["jaix.ExperimentConfig"]["opt_config"] = {
            "jaix.runner.ask_tell.ATOptimiserConfig": {
                "strategy_class": "jaix.runner.ask_tell.strategy.CMA",
                "strategy_config": {
                    "jaix.runner.ask_tell.strategy.CMAConfig": {
                        "sigma0": 2,
                    },
                },
                "init_pop_size": 1,
                "stop_after": 400,
            },
        }
    print(xconfig)
    return xconfig


def test_wandb_logger():
    exp_config = CF.from_dict(get_config())
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
    run = wandb_init(run_config=deepcopy(get_config()), project="ci-cd")
    assert run.mode == "dryrun"
    shutil.rmtree(run.dir)
    run.finish()

    os.environ["WANDB_MODE"] = prev_mode


@pytest.mark.parametrize(
    "suite, comp", itertools.product(["COCO", "RBF", "HPO", "MMind"], [False, True])
)
def test_launch_jaix_experiment(suite, comp):
    prev_mode = os.environ.get("WANDB_MODE", "online")
    os.environ["WANDB_MODE"] = "offline"
    data_dir, exit_code = launch_jaix_experiment(
        run_config=deepcopy(get_config(suite, comp)), project="ci-cd"
    )

    # Remove logging files
    shutil.rmtree(data_dir)
    os.environ["WANDB_MODE"] = prev_mode

    assert exit_code == 0
