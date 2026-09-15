import json
import os
import uuid

import numpy as np
import pandas as pd

from config_postprocess import parse_args
from plots_grid import plot_grid
from plots_sankey import plot_sankey_flows
from plots_scatter import plot_scatter
from utils_read import get_nsga3x_results
from utils_space import angle_between


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
    output_dir: str | None = None,
    file_prefix: str | None = None,
) -> tuple[pd.DataFrame, str | None]:
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
    if output_dir is not None and file_prefix is not None:
        # save dataframe to csv
        file_name = (
            f"{file_prefix}_per_niche_success.csv"
            if per_niche
            else f"{file_prefix}_niche_success.csv"
        )
        file_path = os.path.join(output_dir, file_name)
        niche_success.to_csv(file_path, index=False)
    else:
        file_path = None

    return niche_success, file_path


def compile_distance_success(
    data: pd.DataFrame,
    num_parents: int = 2,
    ideal_point: np.ndarray | None = None,
    output_dir: str | None = None,
    file_prefix: str | None = None,
) -> tuple[pd.DataFrame, str | None]:
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

    # pick the relevant columns for distance success analysis
    distance_success = df[additional_cols + success_cols]
    if output_dir is not None and file_prefix is not None:
        file_name = f"{file_prefix}_distance_success.csv"
        file_path = os.path.join(output_dir, file_name)
        distance_success.to_csv(file_path, index=False)
    else:
        file_path = None
    return distance_success, file_path


def plot_sankey(source_df: pd.DataFrame, output_dir: str, file_prefix: str) -> str:
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
    return plot_sankey_flows(
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
    results_dir: str | None = None,
    problem_name: str | None = None,
) -> tuple[pd.DataFrame, str | None]:
    """
    Compile the data for prediction of offspring success based on parent distances and niches.
    """
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
    if results_dir is not None and problem_name is not None:
        file_name = f"{problem_name}_pred_data.csv"
        file_path = os.path.join(results_dir, file_name)

        pred_data.to_csv(
            os.path.join(results_dir, f"{problem_name}_pred_data.csv"), index=False
        )
    else:
        file_path = None
    return pred_data, file_path


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
    p_dict: dict, out_dir: str = ".", skip_plots: bool = False
) -> dict:

    result_files: dict[int, dict] = {}

    for problem, problem_dict in p_dict.items():
        problem_name = problem_dict.get("problem", problem)
        print(f"Processing results for problem {problem}: {problem_name}")
        result_files[problem] = {}

        runs = [key for key in problem_dict if key.startswith("r_")]
        res_files = [problem_dict[run]["result_file"] for run in runs]

        # combine all the result files into a single dataframe
        problem_df = pd.concat([pd.read_csv(f) for f in res_files], ignore_index=True)

        # Create dataframes for niche success and distance success
        per_niche_success, pniche_csv = compile_niche_success(
            problem_df,
            num_parents=2,
            per_niche=True,
            output_dir=out_dir,
            file_prefix=problem_name,
        )

        niche_success, niche_csv = compile_niche_success(
            problem_df,
            num_parents=2,
            per_niche=False,
            output_dir=out_dir,
            file_prefix=problem_name,
        )

        distance_success, dist_csv = compile_distance_success(
            problem_df,
            num_parents=2,
            ideal_point=problem_dict["ideal_point"],
            output_dir=out_dir,
            file_prefix=problem_name,
        )
        _, pred_csv = compile_pred_data(
            distance_success, problem_df, results_dir=out_dir, problem_name=problem_name
        )
        result_files[problem]["data"] = {
            "per_niche_success": pniche_csv,
            "niche_success": niche_csv,
            "distance_success": dist_csv,
            "pred_data": pred_csv,
        }

        if skip_plots:
            continue

        # Plotting
        plot_files: dict[str, str | list[str]] = {}

        plot_files["sankey"] = plot_sankey(
            per_niche_success,
            output_dir=out_dir,
            file_prefix=problem_name,
        )

        # FIXME: this is a workaround because we know that the number of niches is the same for all runs of a problem, but we should probably check that
        num_niches = problem_dict[runs[0]]["config"]["num_offspring"]
        grid_list = []
        for success_col in grid_success_cols:
            gplot_file = plot_grid(
                niche_success,
                max_grid=num_niches,
                grid_colx="parent_0_niche",
                grid_coly="parent_1_niche",
                hue_col=success_col,
                output_dir=out_dir,
                file_prefix=problem_name,
            )
            grid_list.append(gplot_file)
        plot_files["grid"] = grid_list

        scatter_gen_list = []
        scatter_success_list = []
        for success_col in scatter_success_cols:
            for dist_col in ["parent_x_distance", "parent_y_distance", "parent_angle"]:
                splot_file = plot_scatter(
                    distance_success,
                    output_dir=out_dir,
                    x_col=dist_col,
                    y_col=success_col,
                    hue_col="generation",
                    file_prefix=problem_name,
                )
                scatter_gen_list.append(splot_file)

            splot_file = plot_scatter(
                distance_success,
                output_dir=out_dir,
                x_col="parent_x_distance",
                y_col="parent_y_distance",
                hue_col=success_col,
                file_prefix=problem_name,
            )
            scatter_success_list.append(splot_file)

            splot_file = plot_scatter(
                distance_success,
                output_dir=out_dir,
                x_col="parent_x_distance",
                y_col="parent_angle",
                hue_col=success_col,
                file_prefix=problem_name,
            )
            scatter_success_list.append(splot_file)
        plot_files["scatter_gen"] = scatter_gen_list
        plot_files["scatter_success"] = scatter_success_list
        result_files[problem]["plots"] = plot_files

    return result_files


def run_postprocess(
    results_dir: str,
    out_dir: str,
    problem_ids: list[int] | None = None,
    skip_plots: bool = False,
) -> dict:
    results_dict = get_nsga3x_results(results_dir=results_dir, problem_ids=problem_ids)
    exp_id = uuid.uuid4().hex
    output_dir = os.path.join(out_dir, str(exp_id))
    os.makedirs(output_dir, exist_ok=True)
    start_time = pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")
    res_files = postprocess_results(
        results_dict, out_dir=output_dir, skip_plots=skip_plots
    )
    end_time = pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")
    duration = (pd.Timestamp(end_time) - pd.Timestamp(start_time)).total_seconds()

    for problem_id, problem_dict in results_dict.items():
        problem_dict.update(res_files[problem_id])

    # save results_dict to json file
    file_name = f"results_dict_{exp_id}.json"
    file_path = os.path.join(output_dir, file_name)

    meta_dict = {
        "exp_id": exp_id,
        "start_time": start_time,
        "end_time": end_time,
        "duration_seconds": duration,
        "res_file": file_path,
    }
    results_dict["meta"] = meta_dict
    with open(file_path, "w") as f:
        json.dump(results_dict, f, indent=4, default=str)

    return results_dict


def main(args):
    res = run_postprocess(
        results_dir=args.results_dir,
        out_dir=args.out_dir,
        problem_ids=args.problem_ids,
        skip_plots=args.skip_plots,
    )
    return res


if __name__ == "__main__":
    args = parse_args()
    main(args)
