from jaix import Experiment
from ttex.config import ConfigFactory as CF
from wandb.sdk import launch
from importlib.metadata import version
from typing import Dict
import os
import wandb
from wandb import AlertLevel


def launch_jaix_experiment(run_config: Dict):
    # Update config
    jaix_version = version("tai_jaix")
    config_override = {"repo": "jaix", "version": jaix_version}

    run_config.update(config_override)
    run = wandb.init(config=run_config)

    exp_config = CF.from_dict(run_config)

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
    # TODO: split up prepare and run to add wandbhandler
    # and logging wrapper?


if __name__ == "__main__":
    # This is to test launch from wandb
    if not os.environ.get("WANDB_CONFIG", None):
        raise RuntimeError("Needs to be launched from wandb")
    run_config = launch.load_wandb_config()
    launch_jaix_experiment(run_config)
