import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import (
    HistGradientBoostingClassifier,
    HistGradientBoostingRegressor,
)
from sklearn.feature_selection import mutual_info_classif, mutual_info_regression
from sklearn.inspection import permutation_importance
from sklearn.model_selection import (
    LeaveOneGroupOut,
    cross_val_score,
)

logger = logging.getLogger(__name__)


def get_model(
    max_iter: int = 200,
    learning_rate: float = 0.05,
    max_leaf_nodes: int = 15,
    l2_regularization: float = 1.0,
    random_state: int | None = None,
    target_type: str = "binary",
    **kwargs,
):
    if target_type == "binary":
        model = HistGradientBoostingClassifier(
            max_iter=max_iter,
            learning_rate=learning_rate,
            max_leaf_nodes=max_leaf_nodes,
            l2_regularization=l2_regularization,
            random_state=random_state,
        )
    elif target_type in ("ordinal", "regression"):
        model = HistGradientBoostingRegressor(
            max_iter=max_iter,
            learning_rate=learning_rate,
            max_leaf_nodes=max_leaf_nodes,
            l2_regularization=l2_regularization,
            random_state=random_state,
        )
    else:
        raise ValueError("target_type must be 'binary', 'ordinal', or 'regression'")
    return model


def get_mutual_information(
    X: pd.DataFrame, y: pd.Series, target_type: str, random_state: int | None = None
):
    if target_type == "binary":
        mi = mutual_info_classif(
            X,
            y,
            random_state=random_state,
        )
    elif target_type in ("ordinal", "regression"):
        mi = mutual_info_regression(
            X,
            y,
            random_state=random_state,
        )
    else:
        raise ValueError("target_type must be 'binary', 'ordinal', or 'regression'")
    return mi


def get_scoring_method(
    target_type: str = "binary",
    scoring: tuple[str, bool] | None = None,
) -> tuple[str, bool]:
    if scoring is not None:
        return scoring
    elif target_type == "binary":
        return "roc_auc", True
    elif target_type in ("ordinal", "regression"):
        return "r2", True
    else:
        raise ValueError("target_type must be 'binary', 'ordinal', or 'regression'")


def get_cv_score(
    X: pd.DataFrame,
    y: pd.Series,
    model,
    cv,
    scoring: str,
    groups: pd.Series,
    input_cols: list[str] | None = None,
):
    if input_cols is not None:
        X_to_use = X[input_cols]
    else:
        X_to_use = X

    cv_scores = cross_val_score(
        model,
        X_to_use,
        y,
        cv=cv,
        groups=groups,
        scoring=scoring,
        n_jobs=-2,  # Use all CPUs except one to avoid overloading the system
        error_score="raise",
    )

    return cv_scores.mean(), cv_scores.std()


def permutation_importance_analysis(
    X: pd.DataFrame,
    y: pd.Series,
    model,
    cv,
    scoring: str,
    groups: pd.Series,
    n_permutation_repeats: int = 10,
    random_state: int | None = None,
    **kwargs,
):

    importance_vals: list[np.ndarray] = []
    seeds = np.random.RandomState(random_state).randint(
        0, 10000, size=cv.get_n_splits(X, y, groups=groups)
    )

    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X, y, groups=groups)):

        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]

        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]

        fitted_model = clone(model)
        fitted_model.fit(X_train, y_train)

        perm = permutation_importance(
            fitted_model,
            X_test,
            y_test,
            scoring=scoring,
            n_repeats=n_permutation_repeats,
            random_state=seeds[fold_idx],
        )
        avg_importance = np.mean(perm.importances, axis=1)
        importance_vals.append(avg_importance)

    mean_perm_importance = np.mean(importance_vals, axis=0)
    std_perm_importance = np.std(importance_vals, axis=0)
    return mean_perm_importance, std_perm_importance


def get_loo_cv_scores(
    X: pd.DataFrame,
    y: pd.Series,
    model,
    cv,
    scoring: str,
    input_cols: list[str],
    groups: pd.Series,
) -> tuple[list[float], list[float]]:
    loo_cv_scores = []
    loo_cv_scores_std = []
    for input_col in input_cols:
        # leave-one-feature-out analysis
        input_cols_reduced = [col for col in input_cols if col != input_col]
        loo_cv_score_mean, loo_cv_score_std = get_cv_score(
            X,
            y,
            model,
            cv,
            scoring,
            input_cols=input_cols_reduced,
            groups=groups,
        )
        loo_cv_scores.append(loo_cv_score_mean)
        loo_cv_scores_std.append(loo_cv_score_std)
    return loo_cv_scores, loo_cv_scores_std


@dataclass
class PredResult:
    cv_score_mean: float
    cv_score_std: float
    perm_importance: pd.Series | None
    mutual_info: pd.Series | None
    loo_cv_scores: pd.Series | None
    loo_cv_scores_std: pd.Series | None
    loo_score_drop: pd.Series | None
    loo_score_drop_rel: pd.Series | None


def score_drop(
    values: pd.Series, baseline: float, larger_is_better: bool = True
) -> pd.Series:
    if larger_is_better:
        return pd.Series(baseline - values)
    else:
        return pd.Series(values - baseline)


def run_analysis(
    df: pd.DataFrame,
    input_cols: list[str],
    target_col: str,
    target_type: str,
    group_cols: list[str],
    skip_feature_analysis: bool = False,
    random_state: int | None = None,
    **kwargs,
) -> PredResult:
    X = df[input_cols]
    y = df[target_col]
    logger.info(
        f"Running analysis for target {target_col} with target type {target_type}"
    )
    logger.debug(f"Data shape: {X.shape}, Target shape: {y.shape}")

    model = get_model(target_type=target_type, **kwargs)
    cv = LeaveOneGroupOut()
    scoring, larger_is_better = get_scoring_method(target_type=target_type)
    groups = df[group_cols].astype(str).agg("_".join, axis=1)

    cv_score_mean, cv_score_std = get_cv_score(
        X, y, model, cv, scoring, input_cols=input_cols, groups=groups
    )
    logger.debug(
        f"Cross-validation score for target {target_col}: {cv_score_mean:.4f} ± {cv_score_std:.4f}"
    )
    if skip_feature_analysis:
        logger.info(f"Skipping feature analysis for target {target_col}")
        return PredResult(
            cv_score_mean=cv_score_mean,
            cv_score_std=cv_score_std,
            perm_importance=None,
            mutual_info=None,
            loo_cv_scores=None,
            loo_cv_scores_std=None,
            loo_score_drop=None,
            loo_score_drop_rel=None,
        )

    mutual_info = get_mutual_information(
        X, y, target_type=target_type, random_state=random_state
    )
    logger.debug(f"Mutual information for target {target_col}:")
    logger.debug(mutual_info)

    perm_importance, _perm_importance_std = permutation_importance_analysis(
        X, y, model, cv, scoring, groups=groups, **kwargs
    )
    logger.debug(f"Permutation importance for target {target_col} computed")
    logger.debug(perm_importance)

    loo_cv_scores, loo_cv_scores_std = get_loo_cv_scores(
        X, y, model, cv, scoring, input_cols=input_cols, groups=groups
    )
    logger.debug(
        f"Leave-one-feature-out cross-validation scores for target {target_col}"
    )
    logger.debug(loo_cv_scores)
    res = PredResult(
        cv_score_mean=cv_score_mean,
        cv_score_std=cv_score_std,
        perm_importance=pd.Series(perm_importance, index=input_cols),
        mutual_info=pd.Series(mutual_info, index=input_cols),
        loo_cv_scores=pd.Series(loo_cv_scores, index=input_cols),
        loo_cv_scores_std=pd.Series(loo_cv_scores_std, index=input_cols),
        loo_score_drop=score_drop(
            pd.Series(loo_cv_scores, index=input_cols), cv_score_mean, larger_is_better
        ),
        loo_score_drop_rel=score_drop(
            pd.Series(loo_cv_scores, index=input_cols) / cv_score_mean,
            1.0,
            larger_is_better,
        ),
    )
    return res
