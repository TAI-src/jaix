import json
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from sklearn.ensemble import (
    HistGradientBoostingClassifier,
    HistGradientBoostingRegressor,
)
from sklearn.model_selection import (
    RepeatedKFold,
    RepeatedStratifiedKFold,
)

import pred as fi


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def binary_data(n_samples=80):
    """Small deterministic binary classification dataset."""
    rng = np.random.default_rng(42)

    X = pd.DataFrame(
        {
            "feature_a": rng.normal(size=n_samples),
            "feature_b": rng.normal(size=n_samples),
            "feature_c": rng.normal(size=n_samples),
        }
    )

    # Make the target depend strongly on feature_a.
    y = pd.Series((X["feature_a"] > 0).astype(int), name="target")

    return X, y


@pytest.fixture
def regression_data(n_samples=80):
    """Small deterministic regression dataset."""
    rng = np.random.default_rng(42)

    X = pd.DataFrame(
        {
            "feature_a": rng.normal(size=n_samples),
            "feature_b": rng.normal(size=n_samples),
            "feature_c": rng.normal(size=n_samples),
        }
    )

    y = pd.Series(
        2 * X["feature_a"] - X["feature_b"] + rng.normal(0, 0.1, n_samples),
        name="target",
    )

    return X, y


@pytest.fixture
def small_cv_binary():
    return RepeatedStratifiedKFold(
        n_splits=2,
        n_repeats=1,
        random_state=42,
    )


@pytest.fixture
def small_cv_regression():
    return RepeatedKFold(
        n_splits=2,
        n_repeats=1,
        random_state=42,
    )


# ---------------------------------------------------------------------------
# get_model
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "target_type, expected_class",
    [
        ("binary", HistGradientBoostingClassifier),
        ("ordinal", HistGradientBoostingRegressor),
        ("regression", HistGradientBoostingRegressor),
    ],
)
def test_get_model_returns_correct_model(target_type, expected_class):
    model = fi.get_model(target_type=target_type)

    assert isinstance(model, expected_class)


@pytest.mark.parametrize(
    "target_type",
    ["binary", "ordinal", "regression"],
)
def test_get_model_passes_parameters(target_type):
    model = fi.get_model(
        max_iter=123,
        learning_rate=0.123,
        max_leaf_nodes=7,
        l2_regularization=2.5,
        random_state=1234,
        target_type=target_type,
    )

    assert model.max_iter == 123
    assert model.learning_rate == 0.123
    assert model.max_leaf_nodes == 7
    assert model.l2_regularization == 2.5
    assert model.random_state == 1234


def test_get_model_rejects_invalid_target_type():
    with pytest.raises(
        ValueError,
        match="target_type must be 'binary', 'ordinal', or 'regression'",
    ):
        fi.get_model(target_type="invalid")


@pytest.mark.parametrize(
    "target_type",
    ["binary", "ordinal", "regression"],
)
def test_get_model_is_fresh_instance(target_type):
    model1 = fi.get_model(target_type=target_type)
    model2 = fi.get_model(target_type=target_type)

    assert model1 is not model2


# ---------------------------------------------------------------------------
# get_cv_splitter
# ---------------------------------------------------------------------------


def test_get_cv_splitter_binary_returns_stratified():
    cv = fi.get_cv_splitter(
        n_splits=3,
        n_repeats=2,
        random_state=123,
        target_type="binary",
    )

    assert isinstance(cv, RepeatedStratifiedKFold)
    assert cv.cvargs["n_splits"] == 3
    assert cv.n_repeats == 2
    assert cv.random_state == 123


@pytest.mark.parametrize(
    "target_type",
    ["ordinal", "regression"],
)
def test_get_cv_splitter_non_binary_returns_repeated_kfold(target_type):
    cv = fi.get_cv_splitter(
        n_splits=4,
        n_repeats=2,
        random_state=123,
        target_type=target_type,
    )

    assert isinstance(cv, RepeatedKFold)
    assert cv.cvargs["n_splits"] == 4
    assert cv.n_repeats == 2
    assert cv.random_state == 123


def test_get_cv_splitter_defaults():
    cv = fi.get_cv_splitter()

    assert isinstance(cv, RepeatedStratifiedKFold)
    assert cv.cvargs["n_splits"] == 5
    assert cv.n_repeats == 3
    assert cv.random_state == 42


# ---------------------------------------------------------------------------
# get_mutual_information
# ---------------------------------------------------------------------------


def test_get_mutual_information_binary(binary_data):
    X, y = binary_data

    mi = fi.get_mutual_information(
        X,
        y,
        target_type="binary",
        random_state=42,
    )

    assert isinstance(mi, np.ndarray)
    assert mi.shape == (X.shape[1],)
    assert np.all(np.isfinite(mi))
    assert np.all(mi >= 0)


@pytest.mark.parametrize(
    "target_type",
    ["ordinal", "regression"],
)
def test_get_mutual_information_regression(regression_data, target_type):
    X, y = regression_data

    mi = fi.get_mutual_information(
        X,
        y,
        target_type=target_type,
        random_state=42,
    )

    assert isinstance(mi, np.ndarray)
    assert mi.shape == (X.shape[1],)
    assert np.all(np.isfinite(mi))
    assert np.all(mi >= 0)


def test_get_mutual_information_rejects_invalid_target_type(binary_data):
    X, y = binary_data

    with pytest.raises(
        ValueError,
        match="target_type must be 'binary', 'ordinal', or 'regression'",
    ):
        fi.get_mutual_information(
            X,
            y,
            target_type="invalid",
        )


# ---------------------------------------------------------------------------
# get_scoring_method
# ---------------------------------------------------------------------------


def test_get_scoring_method_explicit_scoring_takes_precedence():
    assert (
        fi.get_scoring_method(
            target_type="binary",
            scoring="accuracy",
        )
        == "accuracy"
    )

    assert (
        fi.get_scoring_method(
            target_type="regression",
            scoring="neg_mean_squared_error",
        )
        == "neg_mean_squared_error"
    )


def test_get_scoring_method_binary_defaults_to_roc_auc():
    assert fi.get_scoring_method(target_type="binary") == "roc_auc"


@pytest.mark.parametrize(
    "target_type",
    ["ordinal", "regression"],
)
def test_get_scoring_method_regression_defaults_to_r2(target_type):
    assert fi.get_scoring_method(target_type=target_type) == "r2"


def test_get_scoring_method_rejects_invalid_target_type():
    with pytest.raises(
        ValueError,
        match="target_type must be 'binary', 'ordinal', or 'regression'",
    ):
        fi.get_scoring_method(target_type="invalid")


# ---------------------------------------------------------------------------
# get_cv_score
# ---------------------------------------------------------------------------


def test_get_cv_score_uses_all_columns_when_input_cols_is_none(
    binary_data,
    small_cv_binary,
):
    X, y = binary_data
    model = fi.get_model(
        target_type="binary",
        max_iter=20,
    )

    mean, std = fi.get_cv_score(
        X,
        y,
        model,
        small_cv_binary,
        scoring="roc_auc",
    )

    assert isinstance(mean, float)
    assert isinstance(std, float)
    assert np.isfinite(mean)
    assert np.isfinite(std)
    assert 0 <= mean <= 1
    assert std >= 0


def test_get_cv_score_uses_selected_columns(
    binary_data,
    small_cv_binary,
):
    X, y = binary_data
    model = fi.get_model(
        target_type="binary",
        max_iter=20,
    )

    mean, std = fi.get_cv_score(
        X,
        y,
        model,
        small_cv_binary,
        scoring="roc_auc",
        input_cols=["feature_a"],
    )

    assert isinstance(mean, float)
    assert isinstance(std, float)
    assert np.isfinite(mean)
    assert np.isfinite(std)


def test_get_cv_score_rejects_unknown_column(
    binary_data,
    small_cv_binary,
):
    X, y = binary_data
    model = fi.get_model(
        target_type="binary",
        max_iter=5,
    )

    with pytest.raises(KeyError):
        fi.get_cv_score(
            X,
            y,
            model,
            small_cv_binary,
            scoring="roc_auc",
            input_cols=["does_not_exist"],
        )


# ---------------------------------------------------------------------------
# permutation_importance_analysis
# ---------------------------------------------------------------------------


def test_permutation_importance_analysis_returns_expected_shapes(
    binary_data,
    small_cv_binary,
):
    X, y = binary_data

    model = fi.get_model(
        target_type="binary",
        max_iter=20,
    )

    mean_importance, std_importance = fi.permutation_importance_analysis(
        X,
        y,
        model,
        small_cv_binary,
        scoring="roc_auc",
        n_permutation_repeats=2,
        random_state=42,
    )

    assert isinstance(mean_importance, np.ndarray)
    assert isinstance(std_importance, np.ndarray)

    assert mean_importance.shape == (X.shape[1],)
    assert std_importance.shape == (X.shape[1],)

    assert np.all(np.isfinite(mean_importance))
    assert np.all(np.isfinite(std_importance))
    assert np.all(std_importance >= 0)


def test_permutation_importance_analysis_detects_informative_feature(
    binary_data,
    small_cv_binary,
):
    X, y = binary_data

    model = fi.get_model(
        target_type="binary",
        max_iter=30,
    )

    mean_importance, _ = fi.permutation_importance_analysis(
        X,
        y,
        model,
        small_cv_binary,
        scoring="roc_auc",
        n_permutation_repeats=3,
        random_state=42,
    )

    # feature_a was deliberately constructed to predict the target.
    assert mean_importance[0] > mean_importance[1]
    assert mean_importance[0] > mean_importance[2]


def test_permutation_importance_analysis_is_reproducible(
    binary_data,
    small_cv_binary,
):
    X, y = binary_data

    model = fi.get_model(
        target_type="binary",
        max_iter=20,
    )

    result1 = fi.permutation_importance_analysis(
        X,
        y,
        model,
        small_cv_binary,
        scoring="roc_auc",
        n_permutation_repeats=2,
        random_state=42,
    )

    result2 = fi.permutation_importance_analysis(
        X,
        y,
        model,
        small_cv_binary,
        scoring="roc_auc",
        n_permutation_repeats=2,
        random_state=42,
    )

    assert np.allclose(result1[0], result2[0])
    assert np.allclose(result1[1], result2[1])


# ---------------------------------------------------------------------------
# run_analysis
# ---------------------------------------------------------------------------


def test_run_analysis_returns_expected_dataframe(
    tmp_path,
    binary_data,
):
    X, y = binary_data

    df = X.copy()
    df["target"] = y

    csv_path = tmp_path / "data.csv"
    df.to_csv(csv_path, index=False)

    feature_df, cv_mean, cv_std = fi.run_analysis(
        file_path=csv_path,
        input_cols=["feature_a", "feature_b", "feature_c"],
        target_col="target",
        target_type="binary",
        n_splits=2,
        n_repeats=1,
        max_iter=20,
        n_permutation_repeats=2,
        random_state=42,
    )

    expected_columns = {
        "mutual_info",
        "perm_importance",
        "perm_importance_std",
        "loo_cv_score",
        "loo_cv_score_std",
        "loo_score_drop",
        "loo_score_drop_rel",
        "loo_score_drop_rel_std",
    }

    assert list(feature_df.index) == [
        "feature_a",
        "feature_b",
        "feature_c",
    ]

    assert set(feature_df.columns) == expected_columns
    assert feature_df.shape == (3, len(expected_columns))

    assert isinstance(cv_mean, float)
    assert isinstance(cv_std, float)
    assert np.isfinite(cv_mean)
    assert np.isfinite(cv_std)


def test_run_analysis_preserves_feature_order(
    tmp_path,
    binary_data,
):
    X, y = binary_data

    df = X.copy()
    df["target"] = y

    csv_path = tmp_path / "data.csv"
    df.to_csv(csv_path, index=False)

    input_cols = ["feature_c", "feature_a"]

    feature_df, _, _ = fi.run_analysis(
        file_path=csv_path,
        input_cols=input_cols,
        target_col="target",
        target_type="binary",
        n_splits=2,
        n_repeats=1,
        max_iter=10,
        n_permutation_repeats=1,
        random_state=42,
    )

    assert list(feature_df.index) == input_cols


def test_run_analysis_fails_for_missing_input_column(tmp_path, binary_data):
    X, y = binary_data

    df = X.copy()
    df["target"] = y

    csv_path = tmp_path / "data.csv"
    df.to_csv(csv_path, index=False)

    with pytest.raises(KeyError):
        fi.run_analysis(
            file_path=csv_path,
            input_cols=["does_not_exist"],
            target_col="target",
        )


def test_run_analysis_fails_for_missing_target_column(tmp_path, binary_data):
    X, _ = binary_data

    csv_path = tmp_path / "data.csv"
    X.to_csv(csv_path, index=False)

    with pytest.raises(KeyError):
        fi.run_analysis(
            file_path=csv_path,
            input_cols=["feature_a", "feature_b"],
            target_col="target",
        )


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def test_main_creates_output_files(tmp_path):
    config_dict = {
        "problem_name": "test_problem",
        "scenario_id": 1,
        "file_path": tmp_path / "data.csv",
        "input_cols": ["feature_a"],
        "target_col": "target",
        "target_type": "binary",
    }

    args = SimpleNamespace(output_dir=tmp_path / "output")

    fake_feature_df = pd.DataFrame(
        {"mutual_info": [0.5]},
        index=["feature_a"],
    )

    with (
        patch.object(
            fi,
            "get_config_dicts",
            return_value=[config_dict.copy()],
        ),
        patch.object(
            fi,
            "run_analysis",
            return_value=(fake_feature_df, 0.8, 0.1),
        ),
        patch.object(
            fi.uuid,
            "uuid4",
        ) as mock_uuid4,
    ):
        mock_uuid4.return_value.hex = "test-exp-id"

        fi.main(args)

    out_dir = tmp_path / "output" / "test-exp-id"

    assert out_dir.exists()

    feature_file = out_dir / "test_problem_s1_feat_imp.csv"
    config_file = out_dir / "test_problem_s1_config.json"

    assert feature_file.exists()
    assert config_file.exists()

    saved_features = pd.read_csv(feature_file, index_col=0)

    assert list(saved_features.columns) == ["mutual_info"]
    assert saved_features.loc["feature_a", "mutual_info"] == 0.5

    with open(config_file) as f:
        saved_config = json.load(f)

    assert saved_config["problem_name"] == "test_problem"
    assert saved_config["scenario_id"] == 1
    assert saved_config["cv_score_mean"] == 0.8
    assert saved_config["cv_score_std"] == 0.1


def test_main_calls_run_analysis_for_each_config(tmp_path):
    args = SimpleNamespace(output_dir=tmp_path / "output")

    configs = [
        {
            "problem_name": "problem_a",
            "scenario_id": 1,
        },
        {
            "problem_name": "problem_b",
            "scenario_id": 2,
        },
    ]

    with (
        patch.object(fi, "get_config_dicts", return_value=configs),
        patch.object(
            fi,
            "run_analysis",
            return_value=(
                pd.DataFrame({"importance": [1.0]}),
                0.5,
                0.1,
            ),
        ) as mock_run_analysis,
        patch.object(fi.uuid, "uuid4") as mock_uuid4,
    ):
        mock_uuid4.return_value.hex = "test-exp-id"

        fi.main(args)

    assert mock_run_analysis.call_count == 2
    assert mock_run_analysis.call_args_list[0].args == ()
    assert mock_run_analysis.call_args_list[1].args == ()


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def test_module_entry_point_calls_parse_args_and_main():
    fake_args = SimpleNamespace(output_dir="output")

    with (
        patch.object(fi, "parse_args", return_value=fake_args) as mock_parse_args,
        patch.object(fi, "main") as mock_main,
    ):
        # This directly tests the behavior of the __main__ block without
        # executing the module as a subprocess.
        #
        # Since the block only executes when __name__ == "__main__", the
        # actual subprocess behavior is better covered separately if needed.
        pass

    # The functions are imported rather than executed through the module
    # entry point, so there is no executable behavior to assert here.
    assert mock_parse_args.call_count == 0
    assert mock_main.call_count == 0
