import argparse
from itertools import product

import numpy as np

from utils_problems import generate_problem_list


def parse_args():

    parser = argparse.ArgumentParser(description="Run NSGA-III on specified problems.")
    parser.add_argument(
        "--n_gen", type=int, required=True, help="Number of generations to run."
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducibility."
    )
    parser.add_argument(
        "--static_ref",
        action="store_true",
        help="Use static reference directions if set.",
    )
    parser.add_argument(
        "--not_static_ref",
        action="store_true",
        help="Use the default dynamic reference directions if set.",
    )
    parser.add_argument(
        "--problem_ids",
        type=int,
        nargs="+",
        default=None,
        help="List of problem IDs to run (0-22).",
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        default=".",
        help="Directory to save the results CSV files.",
    )
    parser.add_argument(
        "--n_runs",
        type=int,
        default=1,
        help="Number of independent runs for each problem.",
    )
    parser.add_argument(
        "--batch_ids",
        type=int,
        nargs="+",
        default=None,
        help="The batch ID for this run to subselect the settings",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed logs during execution.",
    )
    parser.add_argument(
        "--algorithm",
        type=str,
        choices=["nsga2", "nsga3"],
        default="nsga3",
        help="Algorithm to use: 'nsga2' or 'nsga3'. Default is 'nsga3'.",
    )
    args = parser.parse_args()
    if args.static_ref and args.not_static_ref:
        args.static_ref = None  # If both flags are set, treat as None
    elif args.not_static_ref:
        args.static_ref = False
    elif not args.static_ref and not args.not_static_ref:
        args.static_ref = None  # If neither flag is set, treat as None
    return args


def get_batches(
    batch_id: list[int] | None,
    n_runs: int,
    problem_ids: list[int] | None,
    static_ref: bool | None = None,
    seed: int | None = None,
) -> tuple[list[dict], int]:

    problem_list = generate_problem_list(problem_ids)
    rng = np.random.default_rng(seed)
    seeds = rng.integers(0, 1_000_000, size=n_runs)
    static_ref_opts = [True, False] if static_ref is None else [static_ref]
    settings: list[dict] = []
    for problem, run_seed in product(problem_list, seeds):
        for static in static_ref_opts:
            settings.append(
                {
                    "id": len(settings),
                    "seed": run_seed,
                    "problem": problem,
                    "static_ref": static,
                }
            )
    if batch_id is not None:
        return [settings[i] for i in batch_id], len(settings)
    else:
        return settings, len(settings)
