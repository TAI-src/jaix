from jaix import Experiment, LOGGER_NAME
from jaix.experiment import ExperimentConfig
from ttex.config import ConfigFactory as CF
from ttex.log.handler import WandbHandler
from wandb.sdk import launch, AlertLevel
from importlib.metadata import version
from typing import Dict, Optional, List, Any, Tuple
import os
import wandb
from jaix.env.wrapper import LoggingWrapper, LoggingWrapperConfig
from ttex.log import get_logging_config
import sys
import logging
import argparse
import json
from jaix.utils.dict_tools import nested_set
from copy import deepcopy

logger = logging.getLogger(LOGGER_NAME)


def wandb_logger(
    exp_config: ExperimentConfig,
    run: wandb.sdk.wandb_run.Run,
    wandb_logger_name: str = "jaix_wandb",
):
    """
    Add wandb logging to the experiment configuration
    Args:
        exp_config (ExperimentConfig): Experiment configuration
        run (wandb.sdk.wandb_run.Run): Wandb run
        wandb_logger_name (str, optional): Logger name for wandb. Defaults to "jaix_wandb".
    Returns:
        ExperimentConfig: Experiment configuration with wandb logging
    """
    # Adapt LoggingConfig
    if exp_config.logging_config.dict_config:
        logging_config = exp_config.logging_config.dict_config
    else:
        logging_config = get_logging_config(
            logger_name=LOGGER_NAME,
            disable_existing=exp_config.logging_config.disable_existing,
        )
    logging_config["loggers"][wandb_logger_name] = {
        "level": "INFO",
        "handlers": ["wandb_handler"],
    }
    logging_config["handlers"]["wandb_handler"] = {
        "()": WandbHandler,
        "wandb_run": run,
        "custom_metrics": {"env/step": ["env/*"], "restarts/step": ["restarts/*"]},
        "level": "INFO",
    }
    exp_config.logging_config.dict_config = logging_config

    wandb_log_wrapper = (
        LoggingWrapper,
        LoggingWrapperConfig(logger_name=wandb_logger_name),
    )

    if exp_config.env_config.env_wrappers:
        exp_config.env_config.env_wrappers.append(wandb_log_wrapper)
    else:
        exp_config.env_config.env_wrappers = [wandb_log_wrapper]
    return exp_config


def wandb_init(
    run_config: Dict, project: Optional[str] = None, group: Optional[str] = None
):
    """
    Initialize wandb run
    Args:
        run_config (Dict): Run configuration
        project (Optional[str], optional): Wandb project. Defaults to None.
        group (Optional[str], optional): Wandb group. Defaults to None.
    Returns:
        wandb.sdk.wandb_run.Run: Wandb run
    """
    # Config to log
    # jaix_version = version("tai_jaix")
    jaix_version = "0.1.0"  # TODO: Fix this
    config_override = {"repo": "jaix", "version": jaix_version}

    run_config.update(config_override)
    if not project:
        run = wandb.init(config=run_config, group=group)
    else:
        run = wandb.init(config=run_config, project=project, group=group)
    return run


def run_experiment(
    run_config: Dict,
    project: Optional[str] = None,
    wandb: bool = True,
    group_name: Optional[str] = None,
):
    """
    Run an experiment
    Args:
        run_config (Dict): Dictionary with the run configuration
        project (Optional[str], optional): Wandb project. Defaults to None.
        wandb (bool, optional): If True, will log to wandb. Defaults to True.
        group_name (Optional[str], optional): Wandb group name. Defaults to None.
    Returns:
        data_dir (str): Path to the data directory
        exit_code (int): Exit code of the experiment
    """
    run_config = run_config.copy()
    exp_config = CF.from_dict(run_config)
    run = None
    if wandb:
        run = wandb_init(run_config, project=project, group=group_name)
        data_dir = run.dir
        exp_config = wandb_logger(exp_config, run)
        run.alert(
            "Experiment started",
            text="Experiment started",
            level=AlertLevel.INFO,
        )
    else:
        data_dir = None

    try:
        Experiment.run(exp_config)
        exit_code = 0
    except Exception as e:
        logger.error(f"Experiment failed {e}", exc_info=True)
        exit_code = 1

    if run is not None:
        if exit_code == 0:
            run.alert(
                "Experiment ended",
                text="Experiment ended",
                level=AlertLevel.INFO,
            )
        else:
            run.alert(
                "Experiment failed",
                level=AlertLevel.ERROR,
                text="Experiment failed",
            )
        run.finish(exit_code=exit_code)

    return data_dir, exit_code


def launch_jaix_experiment(
    run_config: Dict,
    project: Optional[str] = None,
    wandb: bool = True,
    repeat: int = 1,
    sweep: Optional[Tuple[List[str], List[Any]]] = None,
    group_name: Optional[str] = None,
):
    """
    Launch a jaix experiment from a run_config dictionary
    Args:
        run_config (Dict): Dictionary with the run configuration
        project (Optional[str], optional): Wandb project. Defaults to None.
        wandb (bool, optional): If True, will log to wandb. Defaults to True.
    Returns:
        data_dir (str): Path to the data directory
        exit_code (int): Exit code of the experiment
    """
    run_configs = []
    group_names = []  # type: List[Optional[str]]
    if sweep is not None:
        sweep_keys, sweep_values = sweep
        for sweep_value in sweep_values:
            config = deepcopy(run_config)
            nested_set(config, sweep_keys, sweep_value)
            run_configs.append(deepcopy(config))
            group_names.append(f"{sweep_keys[-1]} {sweep_value}")
    else:
        run_configs.append(run_config)
        group_names.append(None)

    # TODO: make nicer, this just overwrites for now
    if group_name:
        group_names = [group_name] * len(run_configs)
    results = {}

    for run_config, group_name in zip(run_configs, group_names):
        results[group_name] = {
            "run_config": run_config,
            "data_dirs": [],
            "exit_codes": [],
        }
        for _ in range(repeat):
            data_dir, exit_code = run_experiment(run_config, project, wandb, group_name)
            results[group_name]["data_dirs"].append(data_dir)  # type: ignore
            results[group_name]["exit_codes"].append(exit_code)  # type: ignore
    return results


def parse_args():
    parser = argparse.ArgumentParser(description="Launch a jaix experiment")
    parser.add_argument(
        "--project",
        type=str,
        default=None,
        help="Wandb project to log to. If not provided, will not log to wandb",
    )
    parser.add_argument(
        "--config_file", type=str, help="Path to the configuration file"
    )
    parser.add_argument("--repeat", type=int, default=1, help="Number of repetitions")
    parser.add_argument(
        "--sweep_keys", nargs="+", type=str, help="Keys to sweep value in config"
    )
    parser.add_argument("--sweep_values", nargs="+", type=float, help="Values to sweep")
    parser.add_argument("--group", type=str, default=None, help="Wandb group name")
    args = parser.parse_args()
    if args.sweep_values:
        cmp = [int(v) == v for v in args.sweep_values]
        if all(cmp):
            args.sweep_values = [int(v) for v in args.sweep_values]
    return args


if __name__ == "__main__":
    """
    This script is used to launch a jaix experiment from a wandb configuration
    """
    launch_arguments = {}
    if os.environ.get("WANDB_CONFIG", None):
        run_config = launch.load_wandb_config().as_dict()
        launch_arguments["wandb"] = True
        if "repeat" in run_config:
            launch_arguments["repeat"] = run_config.pop("repeat")
        if "sweep" in run_config:
            launch_arguments["sweep"] = run_config.pop("sweep")
        launch_arguments["run_config"] = run_config
    else:
        args = parse_args()
        # run_config = CF.from_file(args.config_file).as_dict()
        with open(args.config_file, "r") as f:
            run_config = json.load(f)
        launch_arguments["run_config"] = run_config
        if args.project:
            launch_arguments["project"] = args.project
            launch_arguments["wandb"] = True
        else:
            launch_arguments["wandb"] = False
        launch_arguments["repeat"] = args.repeat
        if args.sweep_keys and args.sweep_values:
            sweep_keys = args.sweep_keys  # type: List[str]
            sweep_values = args.sweep_values  # type: List[Any]
            launch_arguments["sweep"] = (sweep_keys, sweep_values)  # type: ignore
        launch_arguments["group_name"] = args.group
        # TODO: better validation of arguments
    results = launch_jaix_experiment(**launch_arguments)  # type: ignore
    # Aggregate exit codes. If any experiment failed, the script will return something different than 0
    exit_codes = [max(result["exit_codes"]) for result in results.values()]
    sys.exit(max(exit_codes))
