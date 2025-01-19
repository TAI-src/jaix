import wandb
import numpy as np

api = wandb.Api()
entity, project = "TAI_track", "rbf"
runs = api.runs(entity + "/" + project)

for run in runs:
    """
    importate_opts = run.config["jaix.ExperimentConfig"]["opt_config"][
        "jaix.runner.ask_tell.ATOptimiserConfig"
    ]["strategy_config"]["jaix.runner.ask_tell.strategy.BasicEAConfig"]["update_opts"]
    if "s" not in update_opts:
        update_opts["s"] = np.exp(1) - 1
    run.group = str(update_opts["s"])
    run.update()
    """
    factor = run.config["jaix.ExperimentConfig"]["opt_config"][
        "jaix.runner.ask_tell.ATOptimiserConfig"
    ]["strategy_config"]["jaix.runner.ask_tell.strategy.CMAConfig"]["opts"][
        "popsize_factor"
    ]
    run.group = str(factor)
    run.update()
