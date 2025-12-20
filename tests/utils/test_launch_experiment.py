from jaix.utils.globals import WANDB_LOGGER_NAME
from jaix.utils.launch_experiment import (
    launch_jaix_experiment,
)
import os
from ttex.config import ConfigFactory as CF
from copy import deepcopy
import pytest
import itertools
import json
import logging
from ttex.log import get_wandb_logger
import shutil


def get_config(suite="RBF", comp=False):
    known_suites = {"cont": ["COCO", "RBF", "ASE"], "disc": ["HPO", "MMind"]}

    xconfig = {
        "jaix.experiment.ExperimentConfig": {
            "runner_class": "jaix.runner.ask_tell.ask_tell_runner.ATRunner",
            "runner_config": {
                "jaix.runner.ask_tell.ask_tell_runner.ATRunnerConfig": {
                    "max_evals": 100,
                    "disp_interval": 50,
                },
            },
            "logging_config": {
                "jaix.experiment.LoggingConfig": {
                    "log_level": 10,
                }
            },
        },
    }
    if suite == "COCO":
        xconfig["jaix.experiment.ExperimentConfig"]["env_config"] = {
            "jaix.environment_factory.EnvironmentConfig": {
                "suite_class": "jaix.suite.coco.COCOSuite",
                "suite_config": {
                    "jaix.suite.coco.COCOSuiteConfig": {
                        "env_config": {
                            "jaix.env.singular.ec_env.ECEnvironmentConfig": {
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
        xconfig["jaix.experiment.ExperimentConfig"]["env_config"] = {
            "jaix.environment_factory.EnvironmentConfig": {
                "suite_class": "jaix.suite.ec_suite.ECSuite",
                "suite_config": {
                    "jaix.suite.ec_suite.ECSuiteConfig": {
                        "func_classes": ["jaix.env.utils.problem.rbf_fit.RBFFit"],
                        "func_configs": [
                            {
                                "jaix.env.utils.problem.rbf_fit.RBFFitConfig": {
                                    "rbf_config": {
                                        "jaix.env.utils.problem.rbf.rbf_adapter.RBFAdapterConfig": {},
                                    },
                                    "precision": 1e-3,
                                },
                            }
                        ],
                        "env_config": {
                            "jaix.env.singular.ec_env.ECEnvironmentConfig": {
                                "budget_multiplier": 0.02,
                            },
                        },
                        "instances": list(range(2)),
                        "agg_instances": 1,
                    },
                },
            },
        }
    elif suite == "HPO":
        xconfig["jaix.experiment.ExperimentConfig"]["env_config"] = {
            "jaix.environment_factory.EnvironmentConfig": {
                "suite_class": "jaix.suite.suite.Suite",
                "suite_config": {
                    "jaix.suite.suite.SuiteConfig": {
                        "env_class": "jaix.env.singular.hpo_env.HPOEnvironment",
                        "env_config": {
                            "jaix.env.singular.hpo_env.HPOEnvironmentConfig": {
                                "training_budget": 10,
                                "task_type": "jaix.env.utils.hpo.TaskType.C1",
                                "repo_name": "D244_F3_C1530_30",
                                "cache": True,
                            },
                        },
                        "functions": [0],
                        "agg_instances": 1,
                    },
                },
            },
        }
    elif suite == "MMind":
        xconfig["jaix.experiment.ExperimentConfig"]["env_config"] = {
            "jaix.environment_factory.EnvironmentConfig": {
                "suite_class": "jaix.suite.suite.Suite",
                "suite_config": {
                    "jaix.suite.suite.SuiteConfig": {
                        "env_class": "jaix.env.singular.mastermind_env.MastermindEnvironment",
                        "env_config": {
                            "jaix.env.singular.mastermind_env.MastermindEnvironmentConfig": {
                                "num_slots_range": (4, 6),
                                "num_colours_range": (2, 3),
                                "max_guesses": 10,
                            },
                        },
                        "instances": list(range(2)),
                        "agg_instances": 1,
                    },
                },
            },
        }
    elif suite == "ASE":
        xconfig["jaix.experiment.ExperimentConfig"]["env_config"] = {
            "jaix.environment_factory.EnvironmentConfig": {
                "suite_class": "jaix.suite.suite.Suite",
                "suite_config": {
                    "jaix.suite.suite.SuiteConfig": {
                        "env_class": "jaix.env.singular.LJClustEnvironment",
                        "env_config": {
                            "jaix.env.singular.LJClustEnvironmentConfig": {
                                "ljclust_adapter_config": {
                                    "jaix.env.utils.ase.LJClustAdapterConfig": {
                                        "target_dir": "./tmp_data",
                                    },
                                },
                                "target_accuracy": 0.0,
                            },
                        },
                        "functions": [0],
                        "instances": [0],
                        "agg_instances": 1,
                    },
                },
            },
        }
    xconfig["jaix.experiment.ExperimentConfig"]["env_config"][
        "jaix.environment_factory.EnvironmentConfig"
    ]["seed"] = None
    xconfig["jaix.experiment.ExperimentConfig"]["env_config"][
        "jaix.environment_factory.EnvironmentConfig"
    ]["comp_config"] = None

    if comp:
        xconfig["jaix.experiment.ExperimentConfig"]["env_config"][
            "jaix.environment_factory.EnvironmentConfig"
        ]["comp_config"] = {
            "jaix.environment_factory.CompositeEnvironmentConfig": {
                "agg_type": "jaix.suite.suite.AggType.INST",
                "comp_env_class": "jaix.env.composite.switching_environment.SwitchingEnvironment",
                "comp_env_config": {
                    "jaix.env.composite.switching_environment.SwitchingEnvironmentConfig": {
                        "switching_pattern_class": "jaix.env.utils.switching_pattern.switching_pattern.SeqRegSwitchingPattern",
                        "switching_pattern_config": {
                            "jaix.env.utils.switching_pattern.switching_pattern.SeqRegSwitchingPatternConfig": {
                                "wait_period": 20,
                            },
                        },
                        "real_time": False,
                    },
                },
                "comp_env_wrappers": [
                    (
                        "jaix.env.wrapper.improvement_reward_wrapper.ImprovementRewardWrapper",
                        {
                            "jaix.env.wrapper.improvement_reward_wrapper.ImprovementRewardWrapperConfig": {
                                "state_eval": "obs0"
                            }
                        },
                    )
                ],
            }
        }
        xconfig["jaix.experiment.ExperimentConfig"][
            "opt_class"
        ] = "jaix.runner.ask_tell.at_optimiser.ATOptimiser"
        xconfig["jaix.experiment.ExperimentConfig"]["opt_config"] = {
            "jaix.runner.ask_tell.at_optimiser.ATOptimiserConfig": {
                "strategy_class": "jaix.runner.ask_tell.strategy.bandit.ATBandit",
                "strategy_config": {
                    "jaix.runner.ask_tell.strategy.bandit.ATBanditConfig": {
                        "bandit_config": {
                            "jaix.runner.ask_tell.strategy.utils.bandit_model.BanditConfig": {
                                "epsilon": 0.1,
                                "min_tries": 4,
                                "exploit_strategy": "jaix.runner.ask_tell.strategy.utils.bandit_model.BanditExploitStrategy.MAX",
                            },
                        },
                    },
                },
                "init_pop_size": 1,
            },
        }
        # Continuous optimisation uses CMA-ES, discrete uses Basic EA
        if suite in known_suites["cont"]:
            xconfig["jaix.experiment.ExperimentConfig"]["opt_config"][
                "jaix.runner.ask_tell.at_optimiser.ATOptimiserConfig"
            ]["strategy_config"]["jaix.runner.ask_tell.strategy.bandit.ATBanditConfig"][
                "opt_confs"
            ] = [
                {
                    "jaix.runner.ask_tell.at_optimiser.ATOptimiserConfig": {
                        "strategy_class": "jaix.runner.ask_tell.strategy.cma.CMA",
                        "strategy_config": {
                            "jaix.runner.ask_tell.strategy.cma.CMAConfig": {
                                "sigma0": 2,
                            },
                        },
                        "init_pop_size": 1,
                        "stop_after": 400,
                    }
                },
                {
                    "jaix.runner.ask_tell.at_optimiser.ATOptimiserConfig": {
                        "strategy_class": "jaix.runner.ask_tell.strategy.cma.CMA",
                        "strategy_config": {
                            "jaix.runner.ask_tell.strategy.cma.CMAConfig": {
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
            xconfig["jaix.experiment.ExperimentConfig"]["opt_config"][
                "jaix.runner.ask_tell.at_optimiser.ATOptimiserConfig"
            ]["strategy_config"]["jaix.runner.ask_tell.strategy.bandit.ATBanditConfig"][
                "opt_confs"
            ] = [
                {
                    "jaix.runner.ask_tell.at_optimiser.ATOptimiserConfig": {
                        "strategy_class": "jaix.runner.ask_tell.strategy.basic_ea.BasicEA",
                        "strategy_config": {
                            "jaix.runner.ask_tell.strategy.basic_ea.BasicEAConfig": {
                                "strategy": "jaix.runner.ask_tell.strategy.basic_ea.EAStrategy.Plus",
                                "mu": 1,
                                "lam": 1,
                                "mutation_op": "jaix.runner.ask_tell.strategy.basic_ea.MutationOp.FLIP",
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
                    "jaix.runner.ask_tell.at_optimiser.ATOptimiserConfig": {
                        "strategy_class": "jaix.runner.ask_tell.strategy.basic_ea.BasicEA",
                        "strategy_config": {
                            "jaix.runner.ask_tell.strategy.basic_ea.BasicEAConfig": {
                                "strategy": "jaix.runner.ask_tell.strategy.basic_ea.EAStrategy.Plus",
                                "mu": 2,
                                "lam": 5,
                                "mutation_op": None,
                                "crossover_op": "jaix.runner.ask_tell.strategy.basic_ea.CrossoverOp.UNIFORM",
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
        xconfig["jaix.experiment.ExperimentConfig"][
            "opt_class"
        ] = "jaix.runner.ask_tell.at_optimiser.ATOptimiser"
        if suite in known_suites["cont"]:
            # Continuous optimisation, use CMA-ES
            xconfig["jaix.experiment.ExperimentConfig"]["opt_config"] = {
                "jaix.runner.ask_tell.at_optimiser.ATOptimiserConfig": {
                    "strategy_class": "jaix.runner.ask_tell.strategy.cma.CMA",
                    "strategy_config": {
                        "jaix.runner.ask_tell.strategy.cma.CMAConfig": {
                            "sigma0": 2,
                        },
                    },
                    "init_pop_size": 1,
                    "stop_after": 400,
                },
            }
        else:
            # Discrete optimisation, use BasicEA
            xconfig["jaix.experiment.ExperimentConfig"]["opt_config"] = {
                "jaix.runner.ask_tell.at_optimiser.ATOptimiserConfig": {
                    "strategy_class": "jaix.runner.ask_tell.strategy.basic_ea.BasicEA",
                    "strategy_config": {
                        "jaix.runner.ask_tell.strategy.basic_ea.BasicEAConfig": {
                            "strategy": "jaix.runner.ask_tell.strategy.basic_ea.EAStrategy.Plus",
                            "mu": 1,
                            "lam": 1,
                            "mutation_op": "jaix.runner.ask_tell.strategy.basic_ea.MutationOp.FLIP",
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


def test_integration_wandb_wrapper():
    prev_mode = os.environ.get("WANDB_MODE", "online")
    os.environ["WANDB_MODE"] = "offline"

    xconfig = deepcopy(get_config())
    # Add wandb wrapper
    xconfig["jaix.experiment.ExperimentConfig"]["env_config"][
        "jaix.environment_factory.EnvironmentConfig"
    ]["env_wrappers"] = [
        (
            "jaix.env.wrapper.wandb_wrapper.WandbWrapper",
            {
                "jaix.env.wrapper.wandb_wrapper.WandbWrapperConfig": {
                    "project": "ci-cd",
                    "snapshot": False,
                }
            },
        )
    ]

    results = launch_jaix_experiment(run_config=xconfig)
    exit_code = [result["exit_codes"][0] for result in results.values()][0]

    os.environ["WANDB_MODE"] = prev_mode

    logger = logging.getLogger(WANDB_LOGGER_NAME)
    assert hasattr(logger, "_wandb_setup")
    # This means that the logger was initialised, so everything was activated as planned
    assert get_wandb_logger() is None  # Wandb should be torn down after experiment
    assert exit_code == 0
    shutil.rmtree("./wandb", ignore_errors=True)


def test_integration_coco_wrapper():
    xconfig = deepcopy(get_config())
    # Add coco logger wrapper
    xconfig["jaix.experiment.ExperimentConfig"]["env_config"][
        "jaix.environment_factory.EnvironmentConfig"
    ]["env_wrappers"] = [
        (
            "jaix.env.wrapper.coco_logger_wrapper.COCOLoggerWrapper",
            {"jaix.env.wrapper.coco_logger_wrapper.COCOLoggerWrapperConfig": {}},
        )
    ]

    results = launch_jaix_experiment(run_config=xconfig)
    exit_code = [result["exit_codes"][0] for result in results.values()][0]
    assert exit_code == 0

    data_dir = [result["data_dirs"][0] for result in results.values()][0]
    shutil.rmtree(data_dir, ignore_errors=True)


def check_installed_extras(suite):
    try:
        if suite == "COCO":
            import cocoex  # noqa: F401
        elif suite == "HPO":
            import tabrepo  # noqa: F401
        elif suite == "ASE":
            import ase  # noqa: F401
    except ImportError:
        return False
    return True


@pytest.mark.parametrize(
    "suite, comp",
    itertools.product(["COCO", "RBF", "HPO", "MMind", "ASE"], [False, True]),
)
def test_launch_jaix_experiment(suite, comp):
    config = get_config(suite, comp)
    if not check_installed_extras(suite):
        pytest.skip(
            f"Skipping test for {suite}. Check installed extras if this is unexpected"
        )

    results = launch_jaix_experiment(run_config=deepcopy(config))
    exit_code = [result["exit_codes"][0] for result in results.values()][0]
    assert exit_code == 0


def test_repeat():
    config = get_config("RBF", False)
    results = launch_jaix_experiment(run_config=deepcopy(config), repeat=2)
    exit_codes = [result["exit_codes"] for result in results.values()]
    assert all([len(exit_code) == 2 for exit_code in exit_codes])
    assert all([exit_code == 0 for exit_code in exit_codes[0]])
    assert len(exit_codes) == 1


def test_sweep():
    config = get_config("RBF", False)
    keys = [
        "jaix.experiment.ExperimentConfig",
        "opt_config",
        "jaix.runner.ask_tell.at_optimiser.ATOptimiserConfig",
        "strategy_config",
        "jaix.runner.ask_tell.strategy.cma.CMAConfig",
        "opts",
        "popsize_factor",
    ]
    results = launch_jaix_experiment(
        run_config=deepcopy(config),
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
    if "mmind" in config_file:
        suite = "MMind"
    elif "coco" in config_file:
        suite = "COCO"
    elif "hpo" in config_file:
        suite = "HPO"
    elif "rbf" in config_file:
        suite = "RBF"
    elif "ase" in config_file:
        suite = "ASE"
    else:
        pytest.fail(f"Unknown suite for {config_file}")
    if not check_installed_extras(suite):
        pytest.skip(
            f"Skipping test for {config_file}. Check installed extras if this is unexpected"
        )

    with open(config_file, "r") as f:
        config = json.load(f)
    # modify the config for test (shorter, logging)
    config["jaix.experiment.ExperimentConfig"]["runner_config"][
        "jaix.runner.ask_tell.ask_tell_runner.ATRunnerConfig"
    ]["max_evals"] = 10
    config["jaix.experiment.ExperimentConfig"]["logging_config"][
        "jaix.experiment.LoggingConfig"
    ]["log_level"] = 10
    if (
        config["jaix.experiment.ExperimentConfig"]["env_config"][
            "jaix.environment_factory.EnvironmentConfig"
        ]["suite_class"]
        == "jaix.suite.suite.Suite"
    ):
        config["jaix.experiment.ExperimentConfig"]["env_config"][
            "jaix.environment_factory.EnvironmentConfig"
        ]["suite_config"]["jaix.suite.suite.SuiteConfig"]["functions"] = [0]
        config["jaix.experiment.ExperimentConfig"]["env_config"][
            "jaix.environment_factory.EnvironmentConfig"
        ]["suite_config"]["jaix.suite.suite.SuiteConfig"]["instances"] = [0]
        config["jaix.experiment.ExperimentConfig"]["env_config"][
            "jaix.environment_factory.EnvironmentConfig"
        ]["suite_config"]["jaix.suite.suite.SuiteConfig"]["agg_instances"] = 1

    results = launch_jaix_experiment(run_config=config)

    exit_code = [result["exit_codes"][0] for result in results.values()][0]
    assert exit_code == 0
