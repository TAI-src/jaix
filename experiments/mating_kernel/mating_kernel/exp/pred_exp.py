import argparse
from mating_kernel.exp.exp import Experiment

input_scenarios = {
    "kernel": ["x_dist", "f_dist"],
    "abs": ["p0_X", "p1_x", "p0_F", "p1_F"],
    "pop_fitness": [
        "p1_rank",
        "p1_crowding",
        "p0_rank",
        "p0_crowding",
    ],
    "pop_state": ["b_rank_mean", "b_crowding_mean"],
    # "archive_state": ["b_size", "b_coverage", "b_avg_dist_to_ideal"],
    "age": ["p0_age", "p1_age"],
}


# TODO: compute pangle, pdist, age etc
# Could also use age info at some point
# For now ignore ideal info


class PredExperiment(Experiment):
    @staticmethod
    def parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument(
            "--perf_stats_dir",
            type=str,
            default=None,
            help="Directory to find the performance statistics for the prediction experiment.",
        )
        parser.add_argument(
            "--kernel",
            action="store_true",
            help="Whether to use the kernel for the prediction experiment.",
        )
        parser.add_argument(
            "--use_ideal",
            default=False,
            action="store_true",
            help="Whether to use the ideal point for the prediction experiment.",
        )
        parser.add_argument(
            "--fitness_info",
            action="store_true",
            help="Whether to use the fitness information for the prediction experiment.",
        )
        parser.add_argument(
            "--state",
            action="store_true",
            help="Whether to use the state information for the prediction experiment.",
        )
        parser.add_argument(
            "--target",
            type=str,
            choices=["survived", "o_dist_to_ideal"],
            default="survived",
            help="The target variable for the prediction experiment.",
        )
        return parser
