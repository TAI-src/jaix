import wandb
import numpy as np

api = wandb.Api()
entity, project = "TAI_track", "mmind"
runs = api.runs(entity + "/" + project)

for run in runs:
    update_opts = run.config["jaix.ExperimentConfig"]["opt_config"][
        "jaix.runner.ask_tell.ATOptimiserConfig"
    ]["strategy_config"]["jaix.runner.ask_tell.strategy.BasicEAConfig"]["update_opts"]
    update_opts["group"] = str(update_opts["s"])
    run.group = str(update_opts["s"])
    run.update()
