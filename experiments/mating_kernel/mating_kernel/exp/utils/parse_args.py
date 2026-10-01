import argparse
from mating_kernel.problems.mk_suite import MKSuiteConfig
from mating_kernel.experiments.experiment import ExperimentConfig


def mk_suite_parser():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--cobi", action="store_true", help="Include Cobi problems in the experiment."
    )
    parser.add_argument(
        "--re", action="store_true", help="Include RE problems in the experiment."
    )
    parser.add_argument(
        "--num_objectives",
        type=int,
        nargs="+",
        default=None,
        help="List of number of objectives to include in the experiment.",
    )
    parser.add_argument(
        "--constrained",
        action="store_true",
        help="Include constrained problems in the experiment.",
    )
    return parser


def experiment_parser():
    parser = argparse.ArgumentParser(add_help=False, parents=[mk_suite_parser()])
    parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for the experiment."
    )
    parser.add_argument(
        "--reps", type=int, default=30, help="Number of repetitions for each batch."
    )
    parser.add_argument(
        "--num_batches",
        type=int,
        default=1,
        help="Number of batches to split the experiment into.",
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        required=True,
        help="Output directory for the experiment results.",
    )
    parser.add_argument(
        "--nth", type=int, nargs="+", default=None, help="Run only the nth batch(es)."
    )
    return parser
