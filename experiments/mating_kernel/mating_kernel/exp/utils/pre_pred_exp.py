import pandas as pd
import numpy as np
from pathlib import Path
from mating_kernel.exp.utils.find_files import find_data_files
import logging
import re

logger = logging.getLogger(__name__)


input_scenarios = {
    "kernel": ["p_x_dist", "p_f_dist"],
    "abs": ["p0_X", "p1_X", "p0_F", "p1_F"],
    "rel_fit": [
        "p1_rank",
        "p1_crowding",
        "p0_rank",
        "p0_crowding",
    ],
    "state": ["b_rank_mean", "b_crowding_mean", "b_score"],
    "age": ["p0_age", "p1_age"],
}

target_types = {
    "survived": "binary",
    "o_dist_to_ideal": "regression",
    "o_F_0": "regression",
    "o_F_1": "regression",
}


def add_features(df: pd.DataFrame, feature_names: list[str]) -> pd.DataFrame:
    feature_definitions = {
        "p_x_dist": (
            ["p0_X", "p1_X"],
            lambda df: np.linalg.norm(
                np.stack(df["p0_X"].values) - np.stack(df["p1_X"].values), axis=1
            ),
        ),
        "p_f_dist": (
            ["p0_F", "p1_F"],
            lambda df: np.linalg.norm(
                np.stack(df["p0_F"].values) - np.stack(df["p1_F"].values), axis=1
            ),
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
    missing_features = [
        feature for feature in feature_names if feature not in df.columns
    ]

    if not all(
        feature_name in feature_definitions for feature_name in missing_features
    ):
        logger.info(
            f"One or more unknown features provided: {missing_features}. Known features are: {list(feature_definitions.keys())}"
        )
        # Remove known features from the missing_features list to avoid raising an error for them
        missing_features = [
            feature for feature in missing_features if feature in feature_definitions
        ]

    new_val_dict = {}
    for feature_name in missing_features:
        required_columns, compute = feature_definitions[feature_name]

        missing = [col for col in required_columns if col not in df.columns]
        if missing:
            raise ValueError(
                f"Missing required columns for feature {feature_name}: {missing}"
            )
        new_val_dict[feature_name] = compute(df)
    new_features_df = pd.DataFrame(new_val_dict, index=df.index)

    df = pd.concat([df, new_features_df], axis=1)

    return df


def get_features(
    kernel: bool, abs: bool, rel_fit: bool, state: bool, age: bool
) -> list[str]:
    features = []
    if kernel:
        features.extend(input_scenarios["kernel"])
    if abs:
        features.extend(input_scenarios["abs"])
    if rel_fit:
        features.extend(input_scenarios["rel_fit"])
    if state:
        features.extend(input_scenarios["state"])
    if age:
        features.extend(input_scenarios["age"])
    features.append("seed")  # This is needed done for stratification
    return features


def str_to_np(df: pd.DataFrame, column_names: list[str] | None = None) -> pd.DataFrame:
    if column_names is None:
        # Identify all columns that are of type 'str' and contain array-like strings
        pattern = re.compile(r"^\[.*\]$")
        column_names = [
            col
            for col in df.columns
            if df[col].dtype == "str" and df[col].str.match(pattern).all()
        ]
    new_val_dict = {}
    for column_name in column_names:
        if column_name in df.columns and df[column_name].dtype == "str":
            new_val_dict[column_name] = df[column_name].apply(
                lambda x: np.fromstring(x.strip("[]"), sep=" ")
            )
        else:
            logger.warning(
                f"Column {column_name} is not present in the DataFrame or is not of type 'str', skipping conversion to numpy array."
            )
    if new_val_dict:
        new_features_df = pd.DataFrame(new_val_dict, index=df.index)
        # Remove the original columns that were converted to numpy arrays
        # And concatenate the new features DataFrame with the original DataFrame
        df = pd.concat([df.drop(columns=new_val_dict.keys()), new_features_df], axis=1)
    return df


def expand_array_column(df: pd.DataFrame, column: str) -> pd.DataFrame:
    values = np.stack(df[column].to_numpy())

    expanded = pd.DataFrame(
        values,
        index=df.index,
        columns=[f"{column}_{i}" for i in range(values.shape[1])],
    )
    return expanded


def expand_array_columns(
    df: pd.DataFrame, column_names: list[str] | None = None
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    if column_names is None:
        # Identify all columns that are of type 'object' and contain numpy arrays
        column_names = [
            col
            for col in df.columns
            if df[col].apply(lambda x: isinstance(x, np.ndarray)).all()
        ]
    expanded_dfs = [expand_array_column(df, col) for col in column_names]
    col_names = {col: list(df.columns) for col, df in zip(column_names, expanded_dfs)}
    joined_df = pd.concat([df.drop(columns=column_names)] + expanded_dfs, axis=1)
    return joined_df, col_names


def preprocess(
    df: pd.DataFrame, features: list[str], remove_mutated: bool = True
) -> pd.DataFrame:
    if "seed" not in features:
        logger.warning("Seed column is not in features, adding it for stratification.")
    # remove completely duplicate rows (very unlikely)
    df.drop_duplicates(inplace=True)
    # remove children that were not properly considered (they don't have a generation)
    df.dropna(subset=["o_n_gen"], inplace=True)
    # remove mutated children (i.e. difference between o_X and o_X_b4)
    if remove_mutated:
        df = df[df["o_X"] == df["o_X_b4"]]
    # Turn numpy strings into numpy arrays
    df = str_to_np(df)
    # Add features that are computed from existing columns
    df = add_features(df, features)
    # Expand numpy array columns into separate columns
    df, expanded_names = expand_array_columns(df)
    # Record changes in column names after expansion
    for original_col, new_cols in expanded_names.items():
        if original_col in features:
            features.remove(original_col)
            features.extend(new_cols)

    # Drop all other columns that are not in features
    if not all(feature in df.columns for feature in features):
        missing_features = [
            feature for feature in features if feature not in df.columns
        ]
        raise ValueError(f"Missing required features: {missing_features}")
    df = df[features]
    # Replace inf values with NaN
    df = df.replace([np.inf, -np.inf], np.nan)
    # Drop all rows with NaN values
    df.dropna(inplace=True)
    return df


def get_data(perf_stats_dir: Path | str, pid: str) -> pd.DataFrame:
    perf_merged_files = find_data_files(
        folder=perf_stats_dir,
        file_type_pattern="merged_record_stats_*.csv",  # FIXME: Should not be hardcoded
        problem_names=[pid],
    )
    assert (
        len(perf_merged_files[pid]) == 1
    ), f"Found {len(perf_merged_files[pid])} files for problem {pid}"
    merged_file = perf_merged_files[pid][0]
    df = pd.read_csv(merged_file)
    return df
