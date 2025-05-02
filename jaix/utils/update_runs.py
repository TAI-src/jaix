import wandb
import numpy as np
from collections import defaultdict

api = wandb.Api()
entity, project, group = "TAI_track", "hpo", "all"

if group is not None:
    runs = api.runs(entity + "/" + project, filters={"group": group})
else:
    runs = api.runs(entity + "/" + project)


# agg_instances = defaultdict(int)
max_num = 4
for run in runs:
    print(run.summary)
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
    """
    if "agg_instances" in run.group:
        agg_i = run.config["jaix.ExperimentConfig"]["env_config"][
            "jaix.EnvironmentConfig"
        ]["suite_config"]["jaix.suite.SuiteConfig"]["agg_instances"]
        if agg_instances[agg_i] < max_num:
            run.group = "bandit"
        else:
            run.group = "bandit2"
        agg_instances[agg_i] += 1
        run.update()
    """

# run = runs[0]
# save the metrics for the run to a csv file
# metrics_dataframe = run.history()
# metrics_dataframe.to_csv("metrics.csv")
