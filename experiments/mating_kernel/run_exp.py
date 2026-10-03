import argparse
import logging

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
    parser.add_argument(
        "--log_level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="WARNING",
        help="Set the logging level.",
    )
    args, unknown_args = parser.parse_known_args()

    # Set up logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    experiment_class = EXPERIMENTS[args.experiment]
    experiment_class.run_from_args(unknown_args)


if __name__ == "__main__":
    main()
