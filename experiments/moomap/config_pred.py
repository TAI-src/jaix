import argparse
from itertools import product
from pathlib import Path

import numpy as np

from utils_problems import get_problem_names
from utils_read import find_data_files


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run prediction analysis for MOO-MAP experiments."
    )
    parser.add_argument(
        "--scenario_ids",
        type=int,
        nargs="*",
        help="IDs of the scenario to run (0-31). If not provided, all scenarios will be run.",
        default=None,
        required=False,
    )
    parser.add_argument(
        "--data_dir",
        type=str,
        help="Folder containing the prediction data files.",
        default=str(Path(__file__).parent / "results"),
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Folder to save the feature importance results.",
        default=str(Path(__file__).parent / "pred_results"),
    )
    parser.add_argument(
        "--file_ids",
        type=int,
        nargs="*",
        help="IDs of the data file to run (usually 0-22). If not provided, all files will be run.",
        default=None,
        required=False,
    )
    parser.add_argument(
        "--batch_id",
        type=int,
        nargs="*",
        help="IDs of the batch to run (combination of files and scenarios). If not provided, all combinations will be run.",
        default=None,
        required=False,
    )
    parser.add_argument(
        "--max_iter",
        type=int,
        default=200,
        help="Maximum number of iterations for the model.",
    )
    parser.add_argument(
        "--learning_rate", type=float, default=0.05, help="Learning rate for the model."
    )
    parser.add_argument(
        "--max_leaf_nodes",
        type=int,
        default=15,
        help="Maximum number of leaf nodes for the model.",
    )
    parser.add_argument(
        "--l2_regularization",
        type=float,
        default=1.0,
        help="L2 regularization strength for the model.",
    )
    parser.add_argument(
        "--n_splits", type=int, default=5, help="Number of splits for cross-validation."
    )
    parser.add_argument(
        "--n_repeats",
        type=int,
        default=3,
        help="Number of repeats for cross-validation.",
    )
    parser.add_argument(
        "--n_permutation_repeats",
        type=int,
        default=10,
        help="Number of repeats for permutation importance.",
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducibility."
    )
    args = parser.parse_args()
    return args


def generate_scenario_list():

    target_settings = [
        ("offspring_added", "binary"),
        ("offspring_rank", "ordinal"),
        ("n_offspring_rank", "regression"),
        ("ncrit_offspring_rank", "regression"),
    ]
    archive_stats_cols = [
        "archive_stats_before_coverage",
        "archive_stats_before_unbounded_hv",
        "archive_stats_before_max_rank",
        "archive_stats_before_mean_rank",
    ]
    input_settings = [
        ["parent_0_niche", "parent_1_niche"],
        ["parent_x_distance", "parent_y_distance", "parent_angle"],
        ["parent_x_distance", "parent_0_niche", "parent_1_niche"],
        [
            "parent_x_distance",
            "parent_y_distance",
            "parent_angle",
            "parent_0_dist_to_ideal",
            "parent_1_dist_to_ideal",
        ],
    ]
    scenario_list = []
    for input_cols in input_settings:
        for target_col, target_type in target_settings:
            without_state = {
                "input_cols": input_cols,
                "target_col": target_col,
                "target_type": target_type,
            }
            with_state = {
                "input_cols": input_cols + archive_stats_cols,
                "target_col": target_col,
                "target_type": target_type,
            }
            scenario_list.append(without_state)
            scenario_list.append(with_state)
    return scenario_list


def get_config_dicts(args):
    scenario_list = generate_scenario_list()
    file_list = find_data_files(
        folder=args.data_dir, file_type_pattern="*pred_data.csv"
    )
    scenario_ids = (
        args.scenario_ids
        if args.scenario_ids is not None
        else list(range(len(scenario_list)))
    )
    assert all(0 <= scenario_id < len(scenario_list) for scenario_id in scenario_ids)
    files = [
        (problem_id, file_path)
        for problem_id, file_paths in file_list.items()
        for file_path in file_paths
    ]
    file_ids = args.file_ids if args.file_ids is not None else list(range(len(files)))
    assert all(0 <= file_id < len(files) for file_id in file_ids)
    seed = (
        args.seed
        if args.seed is not None
        else int(np.random.SeedSequence().generate_state(1)[0])
    )

    batch_list = product(scenario_ids, file_ids)
    batch_ids = (
        args.batch_id
        if args.batch_id is not None
        else list(range(len(list(batch_list))))
    )

    config_dicts = []
    for batch_id in batch_ids:
        scenario_id, file_id = list(product(scenario_ids, file_ids))[batch_id]
        scenario = scenario_list[scenario_id]
        problem_id, file_path = files[file_id]
        config_dict = {
            "file_path": file_path,
            "input_cols": scenario["input_cols"],
            "target_col": scenario["target_col"],
            "target_type": scenario["target_type"],
            "max_iter": args.max_iter,
            "learning_rate": args.learning_rate,
            "max_leaf_nodes": args.max_leaf_nodes,
            "l2_regularization": args.l2_regularization,
            "n_splits": args.n_splits,
            "n_repeats": args.n_repeats,
            "n_permutation_repeats": args.n_permutation_repeats,
            "random_state": seed,
            "problem_id": problem_id,
            "problem_name": get_problem_names(problem_ids=[problem_id])[problem_id],
            "file_id": file_id,
            "scenario_id": scenario_id,
            "batch_id": batch_id,
        }
        config_dicts.append(config_dict)
    return config_dicts
