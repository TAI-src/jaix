from jaix import Experiment
from jaix.experiment import ExperimentConfig
from ttex.config import ConfigFactory as CF
from ttex.log.handler import WandbHandler
from wandb.sdk import launch
from importlib.metadata import version
from typing import Dict
import os
import wandb
from wandb import AlertLevel
import logging
from jaix.env.wrapper import LoggingWrapper, LoggingWrapperConfig


def wandb_logger(exp_config: ExperimentConfig, run: wandb.sdk.wandb_run.Run):
    # Set up logger

    logger = logging.getLogger("wandb")
    handler = WandbHandler(run)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)

    wandb_log_wrapper = (LoggingWrapper, LoggingWrapperConfig(logger_name="wandb"))

    if exp_config.env_config.env_wrappers:
        exp_config.env_config.env_wrappers.append(wandb_log_wrapper)
    else:
        exp_config.env_config.env_wrappers = [wandb_log_wrapper]
    return exp_config


def wandb_init(run_config: Dict):
    # Config to log
    jaix_version = version("tai_jaix")
    config_override = {"repo": "jaix", "version": jaix_version}

    run_config.update(config_override)
    run = wandb.init(config=run_config)
    return run


def launch_jaix_experiment(run_config: Dict):
    exp_config = CF.from_dict(run_config)
    run = wandb_init(run_config)
    exp_config = wandb_logger(exp_config, run)

    run.alert("Experiment started", level=AlertLevel.INFO)
    try:
        Experiment.run(exp_config)
    except Exception as e:
        run.alert(
            "Experiment failed",
            level=AlertLevel.ERROR,
            text=str(e),
        )
    run.alert("Experiment ended", level=AlertLevel.INFO)


if __name__ == "__main__":
    # This is to test launch from wandb
    if not os.environ.get("WANDB_CONFIG", None):
        raise RuntimeError("Needs to be launched from wandb")
    run_config = launch.load_wandb_config()
    launch_jaix_experiment(run_config)
