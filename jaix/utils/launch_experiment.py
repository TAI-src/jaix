from jaix.experiment import Experiment
from ttex.config import ConfigFactory as CF
from wandb.sdk import launch
from typing import Dict, Optional, List, Any, Tuple
import os
import sys
import logging
import argparse
import json
from jaix.utils.dict_tools import nested_set
from copy import deepcopy

import jaix.utils.globals as globals

logger = logging.getLogger(globals.LOGGER_NAME)


def run_experiment(
    run_config: Dict,
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
    logger.info(f"Running experiment with config: {exp_config}")
    exp_id = None

    try:
        exp_id = Experiment.run(exp_config)
        logger.info(f"Experiment finished with id: {exp_id}")
        exit_code = 0
    except Exception as e:
        logger.error(f"Experiment failed {e}", exc_info=True)
        exit_code = 1

    return exit_code, exp_id


def launch_jaix_experiment(
    run_config: Dict,
    repeat: int = 1,
    sweep: Optional[Tuple[List[str], List[Any]]] = None,
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
        # If no sweep, just use no group name
        group_names = [None] * len(run_configs)

    results = {}

    for run_config, group_name in zip(run_configs, group_names):
        results[group_name] = {
            "run_config": run_config,
            "data_dirs": [],
            "exit_codes": [],
        }
        # TODO: pass group names through, don't just ignore them
        for _ in range(repeat):
            exit_code, exp_id = run_experiment(run_config)
            results[group_name]["exit_codes"].append(exit_code)  # type: ignore
            results[group_name]["data_dirs"].append(exp_id)  # type: ignore
    return results


def parse_args():
    parser = argparse.ArgumentParser(description="Launch a jaix experiment")
    parser.add_argument(
        "--config_file", type=str, help="Path to the configuration file"
    )
    parser.add_argument("--repeat", type=int, default=1, help="Number of repetitions")
    parser.add_argument(
        "--sweep_keys", nargs="+", type=str, help="Keys to sweep value in config"
    )
    parser.add_argument("--sweep_values", nargs="+", type=float, help="Values to sweep")
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
        launch_arguments["repeat"] = args.repeat
        if args.sweep_keys and args.sweep_values:
            sweep_keys = args.sweep_keys  # type: List[str]
            sweep_values = args.sweep_values  # type: List[Any]
            launch_arguments["sweep"] = (sweep_keys, sweep_values)  # type: ignore
        # TODO: better validation of arguments
    results = launch_jaix_experiment(**launch_arguments)  # type: ignore
    # Aggregate exit codes. If any experiment failed, the script will return something different than 0
    exit_codes = [max(result["exit_codes"]) for result in results.values()]
    sys.exit(max(exit_codes))
