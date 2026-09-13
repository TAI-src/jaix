import json
import os
import uuid
from pathlib import Path

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
    RepeatedKFold,
    RepeatedStratifiedKFold,
    cross_val_score,
)

from config_pred import get_config_dicts, parse_args


def get_model(
    max_iter: int = 200,
    learning_rate: float = 0.05,
    max_leaf_nodes: int = 15,
    l2_regularization: float = 1.0,
    random_state=42,
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


def get_cv_splitter(
    n_splits: int = 5,
    n_repeats: int = 3,
    random_state=42,
    target_type: str = "binary",
    **kwargs,
):
    if target_type == "binary":
        cv = RepeatedStratifiedKFold(
            n_splits=n_splits,
            n_repeats=n_repeats,
            random_state=random_state,
        )
    else:
        cv = RepeatedKFold(
            n_splits=n_splits,
            n_repeats=n_repeats,
            random_state=random_state,
        )
    return cv


def get_mutual_information(
    X: pd.DataFrame, y: pd.Series, target_type: str, random_state=42
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
    scoring: str | None = None,
):
    if scoring is not None:
        return scoring
    elif target_type == "binary":
        return "roc_auc"
    elif target_type in ("ordinal", "regression"):
        return "r2"
    else:
        raise ValueError("target_type must be 'binary', 'ordinal', or 'regression'")


def get_cv_score(
    X: pd.DataFrame,
    y: pd.Series,
    model,
    cv,
    scoring: str,
    input_cols: list[str] | None = None,
):

    if input_cols is not None:
        X_reduced = X[input_cols]

    X_to_use = X_reduced if input_cols is not None else X

    cv_scores = cross_val_score(
        model,
        X_to_use,
        y,
        cv=cv,
        scoring=scoring,
    )

    return cv_scores.mean(), cv_scores.std()


def permutation_importance_analysis(
    X: pd.DataFrame,
    y: pd.Series,
    model,
    cv,
    scoring: str,
    n_permutation_repeats: int = 10,
    random_state: int = 42,
    **kwargs,
):

    importance_vals: list[list[float]] = []

    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X, y)):

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
            random_state=random_state + fold_idx,
        )
        avg_importance = np.mean(perm.importances, axis=1)
        importance_vals.append(avg_importance)

    mean_perm_importance = np.mean(importance_vals, axis=0)
    std_perm_importance = np.std(importance_vals, axis=0)
    return mean_perm_importance, std_perm_importance


def run_analysis(
    file_path: Path,
    input_cols: list[str],
    target_col: str,
    target_type: str = "binary",
    **kwargs,
):
    df = pd.read_csv(file_path)
    X = df[input_cols]
    y = df[target_col]
    print(f"Data shape: {X.shape}, Target shape: {y.shape}")

    model = get_model(target_type=target_type, **kwargs)
    cv = get_cv_splitter(target_type=target_type, **kwargs)
    scoring = get_scoring_method(target_type=target_type)

    feature_df = pd.DataFrame(index=input_cols)
    feature_df["mutual_info"] = get_mutual_information(
        X, y, target_type=target_type, random_state=kwargs.get("random_state", 42)
    )
    print(f"Mutual information for target {target_col}:")
    print(feature_df["mutual_info"])
    feature_df["perm_importance"], feature_df["perm_importance_std"] = (
        permutation_importance_analysis(X, y, model, cv, scoring, **kwargs)
    )
    print(f"Permutation importance for target {target_col}:")
    print(feature_df["perm_importance"])

    cv_score_mean, cv_score_std = get_cv_score(
        X, y, model, cv, scoring, input_cols=input_cols
    )
    print(
        f"Cross-validation score for target {target_col}: {cv_score_mean:.4f} ± {cv_score_std:.4f}"
    )

    loo_cv_scores = []
    loo_cv_scores_std = []
    for input_col in input_cols:
        # leave-one-feature-out analysis
        input_cols_reduced = [col for col in input_cols if col != input_col]
        loo_cv_score_mean, loo_cv_score_std = get_cv_score(
            X, y, model, cv, scoring, input_cols=input_cols_reduced
        )
        loo_cv_scores.append(loo_cv_score_mean)
        loo_cv_scores_std.append(loo_cv_score_std)
    feature_df["loo_cv_score"] = loo_cv_scores
    feature_df["loo_cv_score_std"] = loo_cv_scores_std
    feature_df["loo_score_drop"] = feature_df["loo_cv_score"] - cv_score_mean
    feature_df["loo_score_drop_rel"] = feature_df["loo_score_drop"] / cv_score_mean
    feature_df["loo_score_drop_rel_std"] = feature_df["loo_cv_score_std"] / cv_score_std

    return feature_df, cv_score_mean, cv_score_std


def main(args):

    exp_id = uuid.uuid4().hex
    out_dir = Path(args.output_dir) / exp_id
    os.makedirs(out_dir, exist_ok=True)

    config_dicts = get_config_dicts(args)
    for config_dict in config_dicts:
        print(
            f"Running analysis for problem: {config_dict['problem_name']}, scenario: {config_dict['scenario_id']}"
        )

        print(f"Configuration: {json.dumps(config_dict, indent=4, default=str)}")
        file_prefix = f"{config_dict["problem_name"]}_s{config_dict['scenario_id']}"
        feature_df, cv_score_mean, cv_score_std = run_analysis(**config_dict)
        df_file = f"{file_prefix}_feat_imp.csv"
        feature_df.to_csv(out_dir / df_file)
        print(f"Feature importance saved to {out_dir}/{df_file}")
        config_dict["cv_score_mean"] = cv_score_mean
        config_dict["cv_score_std"] = cv_score_std
        config_file = f"{file_prefix}_config.json"
        with open(out_dir / config_file, "w") as f:
            json.dump(config_dict, f, indent=4, default=str)


if __name__ == "__main__":
    args = parse_args()
    main(args)
