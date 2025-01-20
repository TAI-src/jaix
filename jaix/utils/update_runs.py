import wandb
import numpy as np

api = wandb.Api()
entity, project = "TAI_track", "mmind"
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
    """
    factortor = run.config["jaix.ExperimentConfig"]["opt_config"][
        "jaix.runner.ask_tell.ATOptimiserConfig"
    ]["strategy_config"]["jaix.runner.ask_tell.strategy.CMAConfig"]["opts"][
        "popsize_factor"
    ]
    run.group = str(factor)
    run.update()
    """

run = runs[0]
# save the metrics for the run to a csv file
metrics_dataframe = run.history()
metrics_dataframe.to_csv("metrics.csv")
