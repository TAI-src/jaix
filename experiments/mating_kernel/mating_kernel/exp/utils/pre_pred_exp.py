import pandas as pd
import numpy as np
from pathlib import Path
from mating_kernel.exp.utils.find_files import find_data_files

input_scenarios = {
    "kernel": ["p_x_dist", "p_f_dist"],
    "abs": ["p0_X", "p1_X", "p0_F", "p1_F"],
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

target_types = {
    "survived": "binary",
    "o_dist_to_ideal": "regression",
}


def add_feature(df: pd.DataFrame, feature_name: str) -> pd.DataFrame:
    feature_definitions = {
        "p_x_dist": (
            ["p0_X", "p1_X"],
            lambda df: np.linalg.norm(df["p0_X"].values - df["p1_X"].values, axis=1),
        ),
        "p_f_dist": (
            ["p0_F", "p1_F"],
            lambda df: np.linalg.norm(df["p0_F"].values - df["p1_F"].values, axis=1),
        ),
        "p0_age": (
            ["p0_n_gen", "o_n_gen"],
            lambda df: df["o_n_gen"] - df["p0_n_gen"],
        ),
        "p1_age": (
            ["p1_n_gen", "o_n_gen"],
            lambda df: df["o_n_gen"] - df["p1_n_gen"],
        ),
    }

    if feature_name not in feature_definitions:
        raise ValueError(f"Unknown feature: {feature_name}")

    required_columns, compute = feature_definitions[feature_name]

    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns for feature {feature_name}: {missing}"
        )

    df[feature_name] = compute(df)
    return df


def get_features(
    kernel: bool, abs: bool, pop_fitness: bool, pop_state: bool, age: bool
) -> list[str]:
    features = []
    if kernel:
        features.extend(input_scenarios["kernel"])
    if abs:
        features.extend(input_scenarios["abs"])
    if pop_fitness:
        features.extend(input_scenarios["pop_fitness"])
    if pop_state:
        features.extend(input_scenarios["pop_state"])
    if age:
        features.extend(input_scenarios["age"])
    features.append("seed")  # This is needed done for stratification
    return features


def preprocess(
    df: pd.DataFrame, features: list[str], remove_mutated: bool = True
) -> pd.DataFrame:
    # remove completely duplicate rows (very unlikely)
    df.drop_duplicates(inplace=True)
    # remove children that were not properly considered (they don't have a generation)
    df.dropna(subset=["o_n_gen"], inplace=True)
    # remove mutated children (i.e. difference between o_X and o_X_b4)
    if remove_mutated:
        df = df[df["o_X"] == df["o_X_b4"]]
    for (
        feature
    ) in (
        features
    ):  # Add computed features if they are not already present in the dataframe
        if feature not in df.columns:
            add_feature(df, feature)
    # Drop all other columns that are not in features
    df = df[features]
    return df


def get_data(perf_stats_dir: Path | str, pid: str) -> pd.DataFrame:
    perf_merged_files = find_data_files(
        folder=perf_stats_dir,
        file_type_pattern="merged_record_stats_*_.csv",  # FIXME: Should not be hardcoded
        problem_names=[pid],
    )
    assert (
        len(perf_merged_files) == 1
    ), f"Found {len(perf_merged_files)} files for problem {pid}"
    assert (
        len(perf_merged_files[pid]) == 1
    ), f"Found {len(perf_merged_files[pid])} files for problem {pid}"
    merged_file = perf_merged_files[pid][0]
    df = merged_file.read_csv(merged_file)
    return df
