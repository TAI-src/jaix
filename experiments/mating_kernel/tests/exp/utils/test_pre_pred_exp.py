from mating_kernel.exp.utils.pre_pred_exp import (
    get_data,
    preprocess,
    get_features,
    str_to_np,
    add_features,
    expand_array_column,
    expand_array_columns,
    input_scenarios,
)
from pathlib import Path
import pytest
import pandas as pd
import numpy as np
import logging


@pytest.fixture(scope="module")
def example_data():
    perf_stats_dir = Path(__file__).parent.parent.parent / "data"
    df = get_data(perf_stats_dir, pid="RE21")
    return df


def test_get_data(example_data):
    assert not example_data.empty, "DataFrame should not be empty"


np_cols = ["o_X", "o_F", "p0_X", "p0_F", "p1_X", "p1_F", "o_X_b4"]


@pytest.mark.parametrize("cols", [np_cols, None])
def test_str_to_np(cols, example_data):
    for col in np_cols:
        assert example_data[col].dtype == "str"
    df_converted = str_to_np(example_data.copy(), column_names=cols)
    for col in np_cols:
        assert (
            df_converted[col].dtype == "object"
        ), f"{col} should be of dtype object after conversion"
        assert all(
            isinstance(x, (np.ndarray, np.generic)) for x in df_converted[col]
        ), f"All elements in {col} should be numpy arrays after conversion"


@pytest.mark.parametrize(
    "add_feature_name", ["p_x_dist", "p_f_dist", "p0_age", "p1_age"]
)
def test_add_feature_run(add_feature_name, example_data):
    df = str_to_np(example_data.copy(), column_names=np_cols)
    df_with_feature = add_features(df, [add_feature_name])
    assert (
        add_feature_name in df_with_feature.columns
    ), f"{add_feature_name} should be added to the DataFrame"


def test_add_feature_value():
    df = pd.DataFrame(
        {
            "p0_X": [np.array([0, 0]), np.array([3, 4])],
            "p1_X": [np.array([1, 0]), np.array([7, 8])],
            "p0_F": [np.array([1, 1]), np.array([2, 2])],
            "p1_F": [np.array([1, 3]), np.array([4, 4])],
            "p0_n_gen": [1, 2],
            "o_n_gen": [3, 4],
            "p1_n_gen": [5, 6],
        }
    )
    added_features = ["p_x_dist", "p_f_dist", "p0_age", "p1_age"]
    df_with_features = add_features(df.copy(), list(df.columns) + added_features)

    assert df_with_features["p_x_dist"].iloc[0] == pytest.approx(1)
    assert df_with_features["p_f_dist"].iloc[0] == pytest.approx(2)
    assert df_with_features["p0_age"].iloc[0] == pytest.approx(2)
    assert df_with_features["p1_age"].iloc[0] == pytest.approx(-2)


def test_add_feature_missing_columns(example_data):
    # Remove required columns for p_x_dist
    df_missing_cols = example_data.drop(columns=["p0_X"])
    with pytest.raises(
        ValueError, match="Missing required columns for feature p_x_dist"
    ):
        add_features(df_missing_cols, ["p_x_dist"])


def test_add_feature_unknown_feature(example_data, caplog):
    caplog.set_level(logging.INFO)
    add_features(example_data, ["unknown_feature"])
    assert "One or more unknown features provided" in caplog.text


@pytest.mark.parametrize("remove_mutated", [True, False])
def test_preprocess(remove_mutated, example_data):
    features = list(example_data.columns)  # Use all columns for this test

    preprocessed_df = preprocess(example_data, features)
    assert not preprocessed_df.empty, "Preprocessed DataFrame should not be empty"
    non_np_cols = [col for col in features if col not in np_cols]
    assert all(
        col in preprocessed_df.columns for col in non_np_cols
    ), "All features should be present in the preprocessed DataFrame"
    assert all(f"{col}_0" in preprocessed_df.columns for col in np_cols)

    # Check that there are no NaN values in the required columns
    assert (
        preprocessed_df["o_n_gen"].notna().all()
    ), "There should be no NaN values in the 'o_n_gen' column"
    # Check that there are no infinite values
    assert (
        np.isfinite(preprocessed_df.select_dtypes(include=[np.number])).all().all()
    ), "There should be no infinite values in the DataFrame"


def test_with_computed_features(example_data, caplog):
    features = ["p_x_dist", "p_f_dist", "p0_age", "p1_age"]
    preprocessed_df = preprocess(example_data, features)
    assert all(
        col in preprocessed_df.columns for col in features
    ), "Computed features should be present in the DataFrame"
    assert (
        "seed" not in preprocessed_df.columns
    ), "'seed' should not be present in the DataFrame"
    # Check that there was a warning about the missing seed
    assert caplog.records, "There should be a warning logged about the missing seed"


def test_get_features():
    features = get_features(
        kernel=True, abs=True, rel_fit=True, state=True, age=True, no_x=False
    )
    assert "seed" in features, "'seed' should be included in the features"
    assert isinstance(features, list), "Features should be returned as a list"
    expected_features = [f for flist in input_scenarios.values() for f in flist]
    assert set(features) == set(
        expected_features + ["seed"]
    ), "Features should match expected features"
    features2 = get_features(
        kernel=True, abs=True, rel_fit=True, state=True, age=True, no_x=True
    )
    assert set(features2) == set(features) - {"p_x_dist", "p0_X", "p1_X"}


def test_remove_mutated(caplog):
    caplog.set_level(logging.ERROR)  # Set log high as this is expected to log a warning
    df = pd.DataFrame(
        {
            "o_X": [
                np.array2string(np.array([1, 2])),
                np.array2string(np.array([3, 4])),
                np.array2string(np.array([5, 6])),
                np.array2string(np.array([7, 8])),
            ],
            "o_X_b4": [
                np.array2string(np.array([1, 2])),
                np.array2string(np.array([3, 4])),
                np.array2string(np.array([5, 6.2])),
                np.array2string(np.array([7, 8])),
            ],
            "o_n_gen": [1, 2, 3, 4],
            "a_fevals": [10, 20, 30, 40],
            "o_crowding": [0.1, 0.2, 0.3, 0.4],
            "p1_age": [5, 6, 7, 8],
        }
    )
    preprocessed_df = preprocess(
        df,
        features=["o_X", "o_X_b4", "o_n_gen", "a_fevals", "o_crowding", "p1_age"],
        remove_mutated=True,
    )
    preprocessed_df_no_remove = preprocess(
        df,
        features=["o_X", "o_X_b4", "o_n_gen", "a_fevals", "o_crowding", "p1_age"],
        remove_mutated=False,
    )
    assert len(preprocessed_df) == 3, "One mutated row should be removed"
    assert (
        len(preprocessed_df_no_remove) == 4
    ), "No rows should be removed when remove_mutated is False"


def test_expand_array_column(example_data):
    np_df = str_to_np(example_data.copy(), column_names=["o_X"])
    expanded_df = expand_array_column(np_df, "o_X")
    assert "o_X_0" in expanded_df.columns
    assert "o_X" not in expanded_df.columns
    array_length = len(np_df["o_X"].iloc[0])
    assert len(expanded_df.columns) == array_length


@pytest.mark.parametrize("columns", [["o_X", "p0_X"], None])
def test_expand_array_columns(columns, example_data):
    np_df = str_to_np(example_data.copy(), column_names=columns)
    np_df.dropna(
        subset=["o_n_gen"], inplace=True
    )  # Ensure no NaN values in 'o_n_gen' which generate empty o_F values
    expanded_df, col_names = expand_array_columns(np_df, columns)
    for col in columns or np_cols:
        for i in range(len(np_df[col].iloc[0])):
            assert f"{col}_{i}" in expanded_df.columns
        assert col not in expanded_df.columns
    assert len(example_data.columns) < len(
        expanded_df.columns
    ), "Expanded DataFrame should have more columns than original"

    assert isinstance(col_names, dict)
    assert all(isinstance(v, list) for v in col_names.values())
    assert all(isinstance(k, str) for k in col_names.keys())
    assert set(col_names.keys()) == set(columns or np_cols)
    assert all(
        len(v) == len(np_df[k].iloc[0]) for k, v in col_names.items()
    ), "Each list in col_names should match the length of the corresponding numpy array"
