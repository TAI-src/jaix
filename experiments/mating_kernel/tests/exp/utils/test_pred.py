import numpy as np
import pandas as pd
import pytest
from sklearn.ensemble import (
    HistGradientBoostingClassifier,
    HistGradientBoostingRegressor,
)
from sklearn.model_selection import (
    LeaveOneGroupOut,
)

from mating_kernel.exp.utils.pred import (
    run_analysis,
    get_model,
    get_mutual_information,
    get_scoring_method,
    get_cv_score,
    permutation_importance_analysis,
    PredResult,
    get_loo_cv_scores,
)

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
            "seed": rng.integers(0, 3, size=n_samples),  # 3 groups for LeaveOneGroupOut
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
            "seed": rng.integers(0, 3, size=n_samples),  # 3 groups for LeaveOneGroupOut
        }
    )

    y = pd.Series(
        2 * X["feature_a"] - X["feature_b"] + rng.normal(0, 0.1, n_samples),
        name="target",
    )

    return X, y


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
    model = get_model(target_type=target_type)

    assert isinstance(model, expected_class)


@pytest.mark.parametrize(
    "target_type",
    ["binary", "ordinal", "regression"],
)
def test_get_model_passes_parameters(target_type):
    model = get_model(
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
        get_model(target_type="invalid")


@pytest.mark.parametrize(
    "target_type",
    ["binary", "ordinal", "regression"],
)
def test_get_model_is_fresh_instance(target_type):
    model1 = get_model(target_type=target_type)
    model2 = get_model(target_type=target_type)

    assert model1 is not model2


# ---------------------------------------------------------------------------
# get_mutual_information
# ---------------------------------------------------------------------------


def test_get_mutual_information_binary(binary_data):
    X, y = binary_data

    mi = get_mutual_information(
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

    mi = get_mutual_information(
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
        get_mutual_information(
            X,
            y,
            target_type="invalid",
        )


# ---------------------------------------------------------------------------
# get_scoring_method
# ---------------------------------------------------------------------------


def test_get_scoring_method_explicit_scoring_takes_precedence():
    assert get_scoring_method(
        target_type="binary",
        scoring=("accuracy", True),
    ) == ("accuracy", True)

    assert get_scoring_method(
        target_type="regression",
        scoring=("neg_mean_squared_error", False),
    ) == ("neg_mean_squared_error", False)


def test_get_scoring_method_binary_defaults_to_roc_auc():
    assert get_scoring_method(target_type="binary") == ("roc_auc", True)


@pytest.mark.parametrize(
    "target_type",
    ["ordinal", "regression"],
)
def test_get_scoring_method_regression_defaults_to_r2(target_type):
    assert get_scoring_method(target_type=target_type) == ("r2", True)


def test_get_scoring_method_rejects_invalid_target_type():
    with pytest.raises(
        ValueError,
        match="target_type must be 'binary', 'ordinal', or 'regression'",
    ):
        get_scoring_method(target_type="invalid")


# ---------------------------------------------------------------------------
# get_cv_score
# ---------------------------------------------------------------------------


def test_get_cv_score_uses_all_columns_when_input_cols_is_none(
    binary_data,
):
    X, y = binary_data
    model = get_model(
        target_type="binary",
        max_iter=20,
    )

    mean, std = get_cv_score(
        X, y, model, LeaveOneGroupOut(), scoring="roc_auc", groups=X["seed"]
    )

    assert isinstance(mean, float)
    assert isinstance(std, float)
    assert np.isfinite(mean)
    assert np.isfinite(std)
    assert 0 <= mean <= 1
    assert std >= 0


def test_get_cv_score_uses_selected_columns(
    binary_data,
):
    X, y = binary_data
    model = get_model(
        target_type="binary",
        max_iter=20,
    )

    mean, std = get_cv_score(
        X,
        y,
        model,
        LeaveOneGroupOut(),
        scoring="roc_auc",
        input_cols=["feature_a"],
        groups=X["seed"],
    )

    assert isinstance(mean, float)
    assert isinstance(std, float)
    assert np.isfinite(mean)
    assert np.isfinite(std)


def test_get_cv_score_rejects_unknown_column(
    binary_data,
):
    X, y = binary_data
    model = get_model(
        target_type="binary",
        max_iter=5,
    )

    with pytest.raises(KeyError):
        get_cv_score(
            X,
            y,
            model,
            LeaveOneGroupOut(),
            scoring="roc_auc",
            input_cols=["does_not_exist"],
            groups=X["seed"],
        )


# ---------------------------------------------------------------------------
# permutation_importance_analysis
# ---------------------------------------------------------------------------


def test_permutation_importance_analysis_returns_expected_shapes(
    binary_data,
):
    X, y = binary_data

    model = get_model(
        target_type="binary",
        max_iter=20,
    )

    mean_importance, std_importance = permutation_importance_analysis(
        X,
        y,
        model,
        LeaveOneGroupOut(),
        scoring="roc_auc",
        groups=X["seed"],
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
):
    X, y = binary_data

    model = get_model(
        target_type="binary",
        max_iter=30,
    )

    mean_importance, _ = permutation_importance_analysis(
        X,
        y,
        model,
        LeaveOneGroupOut(),
        scoring="roc_auc",
        n_permutation_repeats=3,
        random_state=42,
        groups=X["seed"],
    )

    # feature_a was deliberately constructed to predict the target.
    assert mean_importance[0] > mean_importance[1]
    assert mean_importance[0] > mean_importance[2]


def test_permutation_importance_analysis_is_reproducible(
    binary_data,
):
    X, y = binary_data

    model = get_model(
        target_type="binary",
        max_iter=20,
    )

    result1 = permutation_importance_analysis(
        X,
        y,
        model,
        LeaveOneGroupOut(),
        scoring="roc_auc",
        n_permutation_repeats=2,
        random_state=42,
        groups=X["seed"],
    )

    result2 = permutation_importance_analysis(
        X,
        y,
        model,
        LeaveOneGroupOut(),
        scoring="roc_auc",
        n_permutation_repeats=2,
        random_state=42,
        groups=X["seed"],
    )

    assert np.allclose(result1[0], result2[0])
    assert np.allclose(result1[1], result2[1])


def test_loo_cv_scores_returns_expected_shape_and_values(
    binary_data,
):
    X, y = binary_data

    model = get_model(target_type="binary")
    cv = LeaveOneGroupOut()

    loo_scores, loo_scores_std = get_loo_cv_scores(
        X,
        y,
        model,
        cv,
        scoring="roc_auc",
        input_cols=["feature_a", "feature_b", "feature_c"],
        groups=X["seed"],
    )
    assert len(loo_scores) == X["seed"].nunique()
    assert len(loo_scores_std) == X["seed"].nunique()
    assert all(isinstance(score, float) for score in loo_scores)


# ---------------------------------------------------------------------------
# run_analysis
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("short", [True, False])
def test_run_analysis_returns_expected_dataframe(
    short,
    binary_data,
):
    X, y = binary_data

    df = X.copy()
    df["target"] = y

    pred_result = run_analysis(
        df,
        input_cols=["feature_a", "feature_b"],
        target_col="target",
        target_type="binary",
        group_cols=["seed"],
        random_state=42,
        skip_feature_analysis=short,
    )
    assert isinstance(pred_result, PredResult)
    assert pred_result.cv_score_mean >= 0
    assert pred_result.cv_score_std >= 0
    if short:
        assert pred_result.perm_importance is None
        assert pred_result.mutual_info is None
        assert pred_result.loo_cv_scores is None
        assert pred_result.loo_cv_scores_std is None
        assert pred_result.loo_score_drop is None
        assert pred_result.loo_score_drop_rel is None
    else:
        assert all(pred_result.loo_score_drop >= 0)
        assert all(pred_result.loo_score_drop_rel >= 0)
        assert all(pred_result.perm_importance >= 0)
        assert all(pred_result.mutual_info >= 0)
        assert len(pred_result.loo_cv_scores) == len(pred_result.loo_cv_scores_std)
        assert len(pred_result.loo_cv_scores) == len(pred_result.loo_score_drop)
        assert len(pred_result.loo_cv_scores) == len(pred_result.loo_score_drop_rel)
