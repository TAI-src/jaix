from jaix import Experiment, LOGGER_NAME
from jaix.experiment import ExperimentConfig
from ttex.config import ConfigFactory as CF
from ttex.log.handler import WandbHandler
from wandb.sdk import launch, AlertLevel
from importlib.metadata import version
from typing import Dict, Optional
import os
import wandb
from jaix.env.wrapper import LoggingWrapper, LoggingWrapperConfig
from ttex.log import get_logging_config
import sys


def wandb_logger(
    exp_config: ExperimentConfig,
    run: wandb.sdk.wandb_run.Run,
    wandb_logger_name: str = "jaix_wandb",
):
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
        "handlers": ["console", "wandb_handler"],
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


def wandb_init(run_config: Dict, project: Optional[str] = None):
    # Config to log
    jaix_version = version("tai_jaix")
    config_override = {"repo": "jaix", "version": jaix_version}

    run_config.update(config_override)
    if not project:
        run = wandb.init(config=run_config)
    else:
        run = wandb.init(config=run_config, project=project)
    return run


def launch_jaix_experiment(run_config: Dict, project: Optional[str] = None):
    exp_config = CF.from_dict(run_config)
    run = wandb_init(run_config, project)
    data_dir = run.dir
    exp_config = wandb_logger(exp_config, run)

    run.alert("Experiment started", text="Experiment started", level=AlertLevel.INFO)
    try:
        Experiment.run(exp_config)
        run.alert("Experiment ended", text="Experiment ended", level=AlertLevel.INFO)
        run.finish(exit_code=0)
    except Exception as e:
        run.alert(
            "Experiment failed",
            level=AlertLevel.ERROR,
            text=str(e),
        )
        run.finish(exit_code=1)
        return data_dir, 1
    return data_dir, 0


if __name__ == "__main__":
    # This is to test launch from wandb
    if not os.environ.get("WANDB_CONFIG", None):
        # TODO: coudl do parseargs in the future
        raise RuntimeError("Needs to be launched from wandb")
    run_config = launch.load_wandb_config().as_dict()
    _, exit_code = launch_jaix_experiment(run_config)
    sys.exit(exit_code)
