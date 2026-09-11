import argparse
import json
import os

import numpy as np
import pandas as pd
from read_utils import build_results_dict
from utils import angle_between

from plots_grid import plot_grid
from plots_sankey import plot_sankey_flows
from plots_scatter import plot_scatter


def get_success_cols(df: pd.DataFrame) -> list:
    """
    Get the columns related to success rate of offspring.
    Adds additional normalised versions of offspring rank as well based on generation
    """
    # normalise offspring rank by generation
    # rank 0 stays 0, worst rank is 1, and all other ranks are scaled between 0 and 1
    df["n_offspring_rank"] = df.groupby(["generation", "seed"])[
        "offspring_rank"
    ].transform(lambda x: x / x.max())
    df.loc[df["offspring_rank"] == 0, "n_offspring_rank"] = 0
    # critical rank is the worst (i.e. max) rank achieved in the generation, which is stored in archive_stats_after_max_rank
    df["critical_rank"] = df["archive_stats_after_max_rank"]

    # normalise offspring rank by critical rank per generation
    df["ncrit_offspring_rank"] = df["offspring_rank"] / df["critical_rank"]
    # if critical rank is 0, set ncrit_offspring_rank to 0
    df.loc[df["critical_rank"] == 0, "ncrit_offspring_rank"] = 0

    success_cols = [
        "offspring_added",
        "offspring_rank",
        "n_offspring_rank",
        "ncrit_offspring_rank",
        "offspring_dist_to_ideal",
    ]
    return success_cols


def compile_niche_success(
    data: pd.DataFrame,
    num_parents: int = 2,
    per_niche: bool = False,
    output_dir: str = ".",
    file_prefix: str = "",
) -> pd.DataFrame:
    """
    Compile the success rate of each niche pair in the data.
    """
    df = data.copy()
    # get different pairings of parent niches and their success rate
    parent_niche_cols = [f"parent_{i}_niche" for i in range(num_parents)]
    # sort the parent niche columns to ensure that the order of parents does not matter
    df[parent_niche_cols] = np.sort(df[parent_niche_cols].values, axis=1)
    if not per_niche:
        group_cols = parent_niche_cols
    else:
        group_cols = parent_niche_cols + ["offspring_niche"]
    success_cols = get_success_cols(df)
    niche_success = (
        df.groupby(group_cols)[success_cols]
        .agg(
            offspring_added_rate=("offspring_added", "mean"),
            offspring_rank_mean=("offspring_rank", "mean"),
            offspring_rank_std=("offspring_rank", "std"),
            num_offspring=("offspring_added", "count"),
            n_offspring_rank_mean=("n_offspring_rank", "mean"),
            n_offspring_rank_std=("n_offspring_rank", "std"),
            ncrit_offspring_rank_mean=("ncrit_offspring_rank", "mean"),
            ncrit_offspring_rank_std=("ncrit_offspring_rank", "std"),
            offspring_dist_to_ideal_mean=(
                "offspring_dist_to_ideal",
                "mean",
            ),  # This auto-removes nans
        )
        .reset_index()
    )

    # save dataframe to csv
    file_name = (
        f"{file_prefix}_per_niche_success.csv"
        if per_niche
        else f"{file_prefix}_niche_success.csv"
    )
    niche_success.to_csv(os.path.join(output_dir, file_name), index=False)

    return niche_success


def compile_distance_success(
    data: pd.DataFrame,
    num_parents: int = 2,
    ideal_point: np.ndarray | None = None,
    output_dir: str = ".",
    file_prefix: str = "",
) -> pd.DataFrame:
    """
    Get data related to success rate of offspring based on the distance between parents.
    """
    assert num_parents == 2, "This function only supports 2 parents."
    df = data.copy()
    # compute parent distances in x and y
    for dim in ["x", "y"]:
        parent_cols = [f"parent_{i}_{dim}" for i in range(num_parents)]
        for col in parent_cols:
            df[col] = df[col].apply(lambda x: np.fromstring(x.strip("[]"), sep=" "))
        df[f"parent_{dim}_distance"] = np.linalg.norm(
            np.stack(df[parent_cols[0]].values) - np.stack(df[parent_cols[1]].values),
            axis=1,
        )

    success_cols = get_success_cols(df)
    additional_cols = [
        "generation",
        "parent_x_distance",
        "parent_y_distance",
        "parent_0_dist_to_ideal",
        "parent_1_dist_to_ideal",
    ]
    # add angle between vectors of parents on y dim with ideal point as origin
    if ideal_point is not None:
        vec1 = np.stack(df["parent_0_y"].values) - ideal_point
        vec2 = np.stack(df["parent_1_y"].values) - ideal_point
        df["parent_angle"] = angle_between(vec1, vec2)
        additional_cols.append("parent_angle")

    # return the relevant columns for distance success analysis
    distance_success = df[additional_cols + success_cols]
    file_name = f"{file_prefix}_distance_success.csv"
    distance_success.to_csv(os.path.join(output_dir, file_name), index=False)
    return distance_success


def plot_sankey(source_df: pd.DataFrame, output_dir: str, file_prefix: str) -> None:
    df = source_df.copy()
    # Prepare the data for the Sankey diagram
    sel_cols = ["parent_0_niche", "parent_1_niche", "offspring_niche", "num_offspring"]
    df = df[sel_cols]
    # remove columns where offspring_niche == -1 (i.e. no offspring added)
    df = df[df["offspring_niche"] != -1]
    # Add parent pair as a new column
    df["parent_pair"] = df.apply(
        lambda row: f"{row['parent_0_niche']}-{row['parent_1_niche']}", axis=1
    )
    plot_sankey_flows(
        df,
        output_dir=output_dir,
        source_col="parent_pair",
        sink_col="offspring_niche",
        value_col="num_offspring",
        file_prefix=file_prefix,
    )


def compile_pred_data(
    distance_success: pd.DataFrame,
    problem_df: pd.DataFrame,
    results_dir: str,
    problem: str,
) -> pd.DataFrame:
    pred_data = distance_success.copy()
    additional_cols = [
        "offspring_niche",
        "parent_0_niche",
        "parent_1_niche",
        "archive_stats_before_coverage",
        "archive_stats_before_unbounded_hv",
        "archive_stats_before_max_rank",
        "archive_stats_before_mean_rank",
        "archive_stats_after_coverage",
        "archive_stats_after_unbounded_hv",
        "archive_stats_after_max_rank",
        "archive_stats_after_mean_rank",
    ]
    pred_data = pred_data.merge(
        problem_df[additional_cols], how="left", left_index=True, right_index=True
    )
    pred_data.to_csv(os.path.join(results_dir, f"{problem}_pred_data.csv"), index=False)
    return pred_data


grid_success_cols = [
    "offspring_added_rate",
    "n_offspring_rank_mean",
    "offspring_dist_to_ideal_mean",
    "num_offspring",
]

scatter_success_cols = [
    "offspring_dist_to_ideal",
    "n_offspring_rank",
    "ncrit_offspring_rank",
]


def postprocess_results(
    p_dict: dict, results_dir: str = ".", skip_plots: bool = False
) -> None:

    for problem, problem_dict in p_dict.items():
        print(f"Processing results for problem: {problem}")

        runs = [key for key in problem_dict if key.startswith("r_")]
        res_files = [problem_dict[run]["file"] for run in runs]
        # FIXME: this is a workaround because we know that the number of niches is the same for all runs of a problem, but we should probably check that
        num_niches = problem_dict[runs[0]]["config"]["num_offspring"]

        # combine all the result files into a single dataframe
        problem_df = pd.concat([pd.read_csv(f) for f in res_files], ignore_index=True)

        # Create data
        per_niche_success = compile_niche_success(
            problem_df,
            num_parents=2,
            per_niche=True,
            output_dir=results_dir,
            file_prefix=problem,
        )

        niche_success = compile_niche_success(
            problem_df,
            num_parents=2,
            per_niche=False,
            output_dir=results_dir,
            file_prefix=problem,
        )

        distance_success = compile_distance_success(
            problem_df, num_parents=2, ideal_point=problem_dict["ideal_point"]
        )
        compile_pred_data(
            distance_success, problem_df, results_dir=results_dir, problem=problem
        )

        if skip_plots:
            continue
        # Plotting

        plot_sankey(
            per_niche_success,
            output_dir=results_dir,
            file_prefix=problem,
        )

        for success_col in grid_success_cols:
            plot_grid(
                niche_success,
                max_grid=num_niches,
                grid_colx="parent_0_niche",
                grid_coly="parent_1_niche",
                success_col=success_col,
                output_dir=results_dir,
                file_prefix=problem,
            )

        for success_col in scatter_success_cols:
            for dist_col in ["parent_x_distance", "parent_y_distance", "parent_angle"]:
                plot_scatter(
                    distance_success,
                    output_dir=results_dir,
                    x_col=dist_col,
                    y_col=success_col,
                    hue_col="generation",
                    file_prefix=problem,
                )
            plot_scatter(
                distance_success,
                output_dir=results_dir,
                x_col="parent_x_distance",
                y_col="parent_y_distance",
                hue_col=success_col,
                file_prefix=problem,
            )

            plot_scatter(
                distance_success,
                output_dir=results_dir,
                x_col="parent_x_distance",
                y_col="parent_angle",
                hue_col=success_col,
                file_prefix=problem,
            )


def parse_args():

    parser = argparse.ArgumentParser(
        description="Postprocess results from NSGA3 experiments."
    )
    parser.add_argument(
        "--results_dir",
        type=str,
        default="results",
        help="Directory containing the results.",
    )
    parser.add_argument(
        "--skip_plots",
        action="store_true",
        help="Skip plotting the results.",
    )
    parser.add_argument(
        "--problem_ids",
        type=int,
        nargs="*",
        default=None,
        help="List of problem IDs to process. If not provided, all problems will be processed.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    results_dict = build_results_dict()
    # save results_dict to json file
    with open("results_dict.json", "w") as f:
        json.dump(results_dict, f, indent=4, default=str)
    postprocess_results(results_dict)
