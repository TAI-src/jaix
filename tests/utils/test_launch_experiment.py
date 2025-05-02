from jaix.utils import launch_jaix_experiment, wandb_logger, wandb_init
import os
import shutil
from ttex.config import ConfigFactory as CF
from copy import deepcopy
import pytest
import itertools
import json


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
            },
        }
    elif suite == "RBF":
        xconfig["jaix.ExperimentConfig"]["env_config"] = {
            "jaix.EnvironmentConfig": {
                "suite_class": "jaix.suite.ECSuite",
                "suite_config": {
                    "jaix.suite.ECSuiteConfig": {
                        "func_classes": ["jaix.env.utils.problem.RBFFit"],
                        "func_configs": [
                            {
                                "jaix.env.utils.problem.RBFFitConfig": {
                                    "rbf_config": {
                                        "jaix.env.utils.problem.rbf.RBFAdapterConfig": {},
                                    },
                                    "precision": 1e-3,
                                },
                            }
                        ],
                        "env_config": {
                            "jaix.env.singular.ECEnvironmentConfig": {
                                "budget_multiplier": 0.02,
                            },
                        },
                        "instances": list(range(2)),
                        "agg_instances": 0,
                    },
                },
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
                                "training_budget": 10,
                                "task_type": "jaix.env.utils.hpo.TaskType.C1",
                                "repo_name": "D244_F3_C1530_30",
                                "cache": True,
                            },
                        },
                        "functions": [0],
                        "agg_instances": 0,
                    },
                },
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
                                "max_guesses": 10,
                            },
                        },
                        "instances": list(range(2)),
                        "agg_instances": 0,
                    },
                },
            },
        }
    xconfig["jaix.ExperimentConfig"]["env_config"]["jaix.EnvironmentConfig"][
        "env_wrappers"
    ] = [("jaix.env.wrapper.AnyFitWrapper", {})]
    xconfig["jaix.ExperimentConfig"]["env_config"]["jaix.EnvironmentConfig"][
        "seed"
    ] = None
    xconfig["jaix.ExperimentConfig"]["env_config"]["jaix.EnvironmentConfig"][
        "comp_config"
    ] = None

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
                                "wait_period": 20,
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
                                "min_tries": 4,
                                "exploit_strategy": "jaix.runner.ask_tell.strategy.utils.BanditExploitStrategy.MAX",
                            },
                        },
                    },
                },
                "init_pop_size": 1,
            },
        }

        if suite == "COCO" or suite == "RBF":
            xconfig["jaix.ExperimentConfig"]["opt_config"][
                "jaix.runner.ask_tell.ATOptimiserConfig"
            ]["strategy_config"]["jaix.runner.ask_tell.strategy.ATBanditConfig"][
                "opt_confs"
            ] = [
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
            ]
        else:
            # Discrete optimisation, use Basic EA
            xconfig["jaix.ExperimentConfig"]["opt_config"][
                "jaix.runner.ask_tell.ATOptimiserConfig"
            ]["strategy_config"]["jaix.runner.ask_tell.strategy.ATBanditConfig"][
                "opt_confs"
            ] = [
                {
                    "jaix.runner.ask_tell.ATOptimiserConfig": {
                        "strategy_class": "jaix.runner.ask_tell.strategy.BasicEA",
                        "strategy_config": {
                            "jaix.runner.ask_tell.strategy.BasicEAConfig": {
                                "strategy": "jaix.runner.ask_tell.strategy.EAStrategy.Plus",
                                "mu": 1,
                                "lam": 1,
                                "mutation_op": "jaix.runner.ask_tell.strategy.MutationOp.FLIP",
                                "crossover_op": None,
                                "mutation_opts": {},
                                "crossover_opts": {},
                            },
                        },
                        "init_pop_size": 1,
                        "stop_after": 400,
                    }
                },
                {
                    "jaix.runner.ask_tell.ATOptimiserConfig": {
                        "strategy_class": "jaix.runner.ask_tell.strategy.BasicEA",
                        "strategy_config": {
                            "jaix.runner.ask_tell.strategy.BasicEAConfig": {
                                "strategy": "jaix.runner.ask_tell.strategy.EAStrategy.Plus",
                                "mu": 2,
                                "lam": 5,
                                "mutation_op": None,
                                "crossover_op": "jaix.runner.ask_tell.strategy.CrossoverOp.UNIFORM",
                                "mutation_opts": {},
                                "crossover_opts": {},
                            },
                        },
                        "init_pop_size": 1,
                        "stop_after": 400,
                    }
                },
            ]

    else:
        xconfig["jaix.ExperimentConfig"][
            "opt_class"
        ] = "jaix.runner.ask_tell.ATOptimiser"
        if suite == "COCO" or suite == "RBF":
            # Continuous optimisation, use CMA-ES
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
        else:
            # Discrete optimisation, use BasicEA
            xconfig["jaix.ExperimentConfig"]["opt_config"] = {
                "jaix.runner.ask_tell.ATOptimiserConfig": {
                    "strategy_class": "jaix.runner.ask_tell.strategy.BasicEA",
                    "strategy_config": {
                        "jaix.runner.ask_tell.strategy.BasicEAConfig": {
                            "strategy": "jaix.runner.ask_tell.strategy.EAStrategy.Plus",
                            "mu": 1,
                            "lam": 1,
                            "mutation_op": "jaix.runner.ask_tell.strategy.MutationOp.FLIP",
                            "crossover_op": None,
                            "mutation_opts": {"p": 0.2},
                            "crossover_opts": {},
                        },
                    },
                    "init_pop_size": 1,
                    "stop_after": 400,
                },
            }
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
    shutil.rmtree(run.dir, ignore_errors=True)
    run.finish()

    os.environ["WANDB_MODE"] = prev_mode


def test_launch_jaix_experiment_wandb():
    prev_mode = os.environ.get("WANDB_MODE", "online")
    os.environ["WANDB_MODE"] = "offline"
    results = launch_jaix_experiment(
        run_config=deepcopy(get_config()), project="ci-cd", wandb=True
    )
    exit_code = [result["exit_codes"][0] for result in results.values()][0]
    data_dir = [result["data_dirs"][0] for result in results.values()][0]
    # Remove logging files
    shutil.rmtree(data_dir, ignore_errors=True)
    os.environ["WANDB_MODE"] = prev_mode

    assert exit_code == 0


@pytest.mark.parametrize(
    "suite, comp", itertools.product(["COCO", "RBF", "HPO", "MMind"], [False, True])
)
def test_launch_jaix_experiment(suite, comp):
    config = get_config(suite, comp)
    results = launch_jaix_experiment(run_config=deepcopy(config), wandb=False)
    exit_code = [result["exit_codes"][0] for result in results.values()][0]
    data_dir = [result["data_dirs"][0] for result in results.values()][0]
    assert data_dir is None
    assert exit_code == 0


def test_repeat():
    config = get_config("COCO", False)
    results = launch_jaix_experiment(run_config=deepcopy(config), wandb=False, repeat=2)
    exit_codes = [result["exit_codes"] for result in results.values()]
    data_dirs = [result["data_dirs"] for result in results.values()]
    assert all([len(exit_code) == 2 for exit_code in exit_codes])
    assert all([len(data_dir) == 2 for data_dir in data_dirs])
    assert all([exit_code == 0 for exit_code in exit_codes[0]])
    assert len(exit_codes) == 1


def test_sweep():
    config = get_config("COCO", False)
    keys = [
        "jaix.ExperimentConfig",
        "opt_config",
        "jaix.runner.ask_tell.ATOptimiserConfig",
        "strategy_config",
        "jaix.runner.ask_tell.strategy.CMAConfig",
        "opts",
        "popsize_factor",
    ]
    results = launch_jaix_experiment(
        run_config=deepcopy(config),
        wandb=False,
        sweep=(keys, [1, 2]),
    )
    exit_codes = [result["exit_codes"][0] for result in results.values()]
    assert all([exit_code == 0 for exit_code in exit_codes])
    assert len(exit_codes) == 2
    assert list(results.keys()) == ["popsize_factor 1", "popsize_factor 2"]


@pytest.mark.parametrize(
    "config_file",
    [
        "/experiments/rbf/brachy.json",
        "/experiments/coco/single_default.json",
        "/experiments/mmind/mmind.json",
        "/experiments/mmind/telltale.json",
        "/experiments/hpo/binary.json",
        "/experiments/mmind/mmind_comp.json",
        "/experiments/hpo/binary_comp.json",
        "/experiments/mmind/telltale_comp.json",
        "/experiments/hpo/binary_warm.json",
        "/experiments/hpo/binary_all.json",
    ],
)
def test_launch_final(config_file):
    with open(config_file, "r") as f:
        config = json.load(f)
    # modify the config for test (shorter, logging)
    config["jaix.ExperimentConfig"]["runner_config"][
        "jaix.runner.ask_tell.ATRunnerConfig"
    ]["max_evals"] = 10
    config["jaix.ExperimentConfig"]["logging_config"]["jaix.LoggingConfig"][
        "log_level"
    ] = 10
    if (
        config["jaix.ExperimentConfig"]["env_config"]["jaix.EnvironmentConfig"][
            "suite_class"
        ]
        == "jaix.suite.Suite"
    ):
        config["jaix.ExperimentConfig"]["env_config"]["jaix.EnvironmentConfig"][
            "suite_config"
        ]["jaix.suite.SuiteConfig"]["functions"] = [0]
        config["jaix.ExperimentConfig"]["env_config"]["jaix.EnvironmentConfig"][
            "suite_config"
        ]["jaix.suite.SuiteConfig"]["instances"] = [0]
        config["jaix.ExperimentConfig"]["env_config"]["jaix.EnvironmentConfig"][
            "suite_config"
        ]["jaix.suite.SuiteConfig"]["agg_instances"] = 0

    results = launch_jaix_experiment(run_config=config, wandb=False)
    exit_code = [result["exit_codes"][0] for result in results.values()][0]
    data_dir = [result["data_dirs"][0] for result in results.values()][0]
    assert data_dir is None
    assert exit_code == 0
