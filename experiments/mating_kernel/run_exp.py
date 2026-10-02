import argparse

from mating_kernel.exp.exp import Experiment
from mating_kernel.exp.perf_exp import PerfExperiment

EXPERIMENTS: dict[str, type[Experiment]] = {
    "perf": PerfExperiment,
}


def main():
    parser = argparse.ArgumentParser(description="Run an experiment.")
    parser.add_argument(
        "experiment", type=str, choices=EXPERIMENTS, help="Experiment to run."
    )
    args, unknown_args = parser.parse_known_args()

    experiment_class = EXPERIMENTS[args.experiment]
    experiment_class.run_from_args(unknown_args)
