import os
from pathlib import Path
from collections import defaultdict
from nsga3_experiment import NSGA3ExperimentConfig
import json
import csv
import pandas as pd
import numpy as np

logs_dir = f"{os.path.dirname(os.path.abspath(__file__))}/logs"
results_dir = f"{os.path.dirname(os.path.abspath(__file__))}/results"

# Get all csv files in the logs directory
csv_files = list(Path(results_dir).rglob("*.csv"))

res_dict = defaultdict(list)
for file in csv_files:
    res_dict[file.name].append(file)

expected_problem_list = NSGA3ExperimentConfig.generate_problem_list()[0:1]


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
    # compute critical rank per generation which is the highest rank per generation where offspring_added is true
    df["critical_rank"] = df.groupby("generation")["offspring_rank"].transform(
        lambda x: x[df.loc[x.index, "offspring_added"]].max()
    )
    # normalise offspring rank by critical rank per generation
    df["ncrit_offspring_rank"] = df["offspring_rank"] / df["critical_rank"]

    success_cols = [
        "offspring_added",
        "offspring_rank",
        "n_offspring_rank",
        "ncrit_offspring_rank",
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
        )
        .reset_index()
    )
    return niche_success


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
        "offspring_dist_to_ideal",
    ]

    # return the relevant columns for distance success analysis
    distance_success = df[additional_cols + success_cols]
    return distance_success


for i, problem in enumerate(expected_problem_list):
    print(
        "Post-processing results for problem: {} with index {}/{}".format(
            problem, i, len(expected_problem_list)
        )
    )
    csv_name = f"results_{problem}.csv"
    res_files = res_dict[csv_name]
    compile_niche_success(pd.read_csv(res_files[0]), num_parents=2)
    compile_distance_success(pd.read_csv(res_files[0]), num_parents=2)
