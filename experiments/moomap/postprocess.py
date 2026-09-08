import os
from pathlib import Path
from collections import defaultdict
from nsga3_experiment import NSGA3ExperimentConfig
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

import json


def build_results_dict(results_dir: str | None = None) -> dict:
    """
    Build a dictionary of results from the logs directory.
    """
    if results_dir is None:
        results_dir = f"{os.path.dirname(os.path.abspath(__file__))}/results"
    res_dict = defaultdict(list)
    # Get all csv files in the logs directory
    csv_files = list(Path(results_dir).rglob("*.csv"))
    for file in csv_files:
        res_dict[file.name].append(file)
    return res_dict


def get_success_cols(df: pd.DataFrame) -> list:
    """
    Get the columns related to success rate of offspring.
    Adds additional normalised versions of offspring rank as well based on generation
    """
    # normalise offspring rank by generation
    # rank 0 stays 0, worst rank is 1, and all other ranks are scaled between 0 and 1
    df["n_offspring_rank"] = df.groupby("generation")["offspring_rank"].transform(
        lambda x: x / x.max()
    )
    df.loc[df["offspring_rank"] == 0, "n_offspring_rank"] = 0
    # compute critical rank per generation which is the highest rank per generation where offspring_added is true
    df["critical_rank"] = df.groupby("generation")["offspring_rank"].transform(
        lambda x: x[df.loc[x.index, "offspring_added"]].max()
    )
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


def compile_niche_success(data: pd.DataFrame, num_parents: int = 2) -> pd.DataFrame:
    """
    Compile the success rate of each niche pair in the data.
    """
    df = data.copy()
    # get different pairings of parent niches and their success rate
    parent_niche_cols = ["parent_{i}_niche".format(i=i) for i in range(num_parents)]
    # sort the parent niche columns to ensure that the order of parents does not matter
    df[parent_niche_cols] = np.sort(df[parent_niche_cols].values, axis=1)
    success_cols = get_success_cols(df)
    niche_success = (
        df.groupby(parent_niche_cols)[success_cols]
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

    return niche_success


def plot_niche_success(
    niche_success: pd.DataFrame,
    max_niche: int,
    output_dir: str,
    success_col: str = "offspring_added_rate",
    file_prefix: str = "",
) -> None:
    """
    Plot the success rate of each niche pair in the data.
    """

    # create a heatmap of the offspring_added_rate for each niche pair
    # niches should be from 0 to max_niche, so we can create a pivot table with all combinations of niches
    idx = pd.MultiIndex.from_product(
        [range(max_niche), range(max_niche)], names=["parent_0_niche", "parent_1_niche"]
    )
    niche_success = (
        niche_success.set_index(["parent_0_niche", "parent_1_niche"])
        .reindex(idx)
        .reset_index()
    )
    pivot_table = niche_success.pivot(
        index="parent_0_niche", columns="parent_1_niche", values=success_col
    )
    cmap = sns.color_palette("viridis", as_cmap=True)
    cmap.set_bad("lightgray")
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        pivot_table,
        annot=False,
        fmt=".2f",
        cmap=cmap,
        cbar_kws={"label": success_col},
    )
    plt.gca().invert_yaxis()
    plt.title(f"{file_prefix}: {success_col} by Parent Niche Pair")
    plt.xlabel("Parent 1 Niche")
    plt.ylabel("Parent 0 Niche")
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, f"{file_prefix}niche_success_{success_col}.png")
    )
    plt.close()


def compile_distance_success(data: pd.DataFrame, num_parents: int = 2) -> pd.DataFrame:
    """
    Get data related to success rate of offspring based on the distance between parents.
    """
    assert num_parents == 2, "This function only supports 2 parents."
    df = data.copy()
    # compute parent distances in x and y
    for dim in ["x", "y"]:
        parent_cols = [
            "parent_{i}_{dim}".format(i=i, dim=dim) for i in range(num_parents)
        ]
        for col in parent_cols:
            df[col] = df[col].apply(lambda x: np.fromstring(x.strip("[]"), sep=" "))
        df["parent_{dim}_distance".format(dim=dim)] = np.linalg.norm(
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

    # return the relevant columns for distance success analysis
    distance_success = df[additional_cols + success_cols]
    return distance_success


def plot_distance_success_2d(
    distance_success: pd.DataFrame,
    output_dir: str,
    success_col: str = "offspring_dist_to_ideal",
    file_prefix: str = "",
) -> None:
    """
    Plot the success rate of offspring based on the distance between parents.
    """
    # filter rows with nan values in success_col
    # This happens because some scores are not computed for offspring if they are not added

    dist_succ = distance_success.copy()
    dist_succ = dist_succ.dropna(subset=[success_col])

    norm = mpl.colors.Normalize(
        vmin=min(dist_succ[success_col]), vmax=max(dist_succ[success_col])
    )
    cmap = plt.get_cmap("viridis")

    plt.figure(figsize=(10, 8))
    ax = sns.scatterplot(
        data=dist_succ,
        x="parent_x_distance",
        y="parent_y_distance",
        hue=success_col,
        palette=cmap,
        hue_norm=norm,
        alpha=0.7,
        legend=False,
    )
    # Create continuous colorbar
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    plt.colorbar(sm, ax=ax, label=success_col)

    only_added = len(distance_success) > len(dist_succ)
    if only_added:
        plt.title(
            f"{file_prefix}: {success_col} by Parent Distance (only offspring added)"
        )
    else:
        plt.title(f"{file_prefix}: {success_col} by Parent Distance")
    plt.xlabel("Parent X Distance")
    plt.ylabel("Parent Y Distance")
    # plt.colorbar(label=success_col)
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, f"{file_prefix}distance_success_{success_col}.png")
    )
    plt.close()


def plot_distance_success_1d(
    distance_success: pd.DataFrame,
    output_dir: str,
    dist_col: str = "parent_x_distance",
    success_col: str = "offspring_dist_to_ideal",
    file_prefix: str = "",
) -> None:
    """
    Plot the success rate of offspring based on the distance between parents in 1D.
    Use the generation as the colour of the points to see if there is a trend over time.
    """
    norm = mpl.colors.Normalize(
        vmin=min(distance_success["generation"]),
        vmax=max(distance_success["generation"]),
    )
    cmap = plt.get_cmap("viridis")
    plt.figure(figsize=(10, 8))
    ax = sns.scatterplot(
        data=distance_success,
        x=dist_col,
        y=success_col,
        hue="generation",
        palette=cmap,
        hue_norm=norm,
        alpha=0.7,
        legend=False,
    )

    # Create continuous colorbar
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    plt.colorbar(sm, ax=ax, label="Generation")

    plt.title(f"{file_prefix}: {success_col} by Parent Distance ({dist_col})")
    plt.xlabel(f"Parent Distance ({dist_col})")
    plt.ylabel(success_col)
    plt.tight_layout()
    plt.savefig(
        os.path.join(
            output_dir,
            f"{file_prefix}_gen_distance_success_{dist_col}_{success_col}.png",
        )
    )
    plt.close()


def postprocess_results(problem_dict: dict, results_dir: str = ".") -> None:

    for problem, res_files in problem_dict.items():
        problem = problem.replace("results_", "").replace(".csv", "")
        print(f"Processing results for problem: {problem}")
        # get example config to know number of niches (dirty)
        # remove results_ prefix and .csv suffix from problem name to get config name
        example_config = res_files[0].parent / f"config_{problem}.json"
        # load the config to get the number of niches

        with open(example_config, "r") as f:
            config = json.load(f)
            num_niches = config["NSGA3ExperimentConfig"]["num_offspring"]

        # combine all the result files into a single dataframe
        problem_df = pd.concat([pd.read_csv(f) for f in res_files], ignore_index=True)
        niche_success = compile_niche_success(problem_df, num_parents=2)
        success_cols = [
            "offspring_added_rate",
            "n_offspring_rank_mean",
            "offspring_dist_to_ideal_mean",
            "num_offspring",
        ]
        for success_col in success_cols:
            plot_niche_success(
                niche_success,
                max_niche=num_niches,
                output_dir=results_dir,
                file_prefix=problem,
                success_col=success_col,
            )

        distance_success = compile_distance_success(problem_df, num_parents=2)
        success_cols = [
            "offspring_dist_to_ideal",
            "n_offspring_rank",
            "ncrit_offspring_rank",
        ]
        for success_col in success_cols:
            for dist_col in ["parent_x_distance", "parent_y_distance"]:
                plot_distance_success_1d(
                    distance_success,
                    output_dir=results_dir,
                    dist_col=dist_col,
                    success_col=success_col,
                    file_prefix=problem,
                )
            plot_distance_success_2d(
                distance_success,
                output_dir=results_dir,
                success_col=success_col,
                file_prefix=problem,
            )


if __name__ == "__main__":
    results_dict = build_results_dict()
    postprocess_results(results_dict)
