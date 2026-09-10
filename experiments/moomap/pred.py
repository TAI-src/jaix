import argparse
import json
import os
import uuid
from itertools import product
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
        importance_vals.append(perm.importances)

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


def generate_scenario_list():

    target_settings = [
        ("offspring_added", "binary"),
        ("offspring_rank", "ordinal"),
        ("n_offspring_rank", "regression"),
        ("n_offspring_rank", "regression"),
    ]
    archive_stats_cols = [
        "archive_stats_before_coverage",
        "archive_stats_before_unbounded_hv",
        "archive_stats_before_max_rank",
        "archive_stats_before_mean_rank",
    ]
    input_settings = [
        ["parent_0_niche", "parent_1_niche"],
        ["parent_x_distance", "parent_y_distance", "parent_angle"],
        ["parent_x_distance", "parent_0_niche", "parent_1_niche"],
        [
            "parent_x_distance",
            "parent_y_distance",
            "parent_angle",
            "parent_0_dist_to_ideal",
            "parent_1_dist_to_ideal",
        ],
    ]
    scenario_list = []
    for input_cols in input_settings:
        for target_col, target_type in target_settings:
            without_state = {
                "input_cols": input_cols,
                "target_col": target_col,
                "target_type": target_type,
            }
            with_state = {
                "input_cols": input_cols + archive_stats_cols,
                "target_col": target_col,
                "target_type": target_type,
            }
            scenario_list.append(without_state)
            scenario_list.append(with_state)
    return scenario_list


def find_data_files(
    folder: str, glob_pattern: str = "*_pred_data.csv"
) -> dict[str, Path]:
    path = Path(folder)
    files = list(path.glob(glob_pattern))
    # get problem names from file names
    problem_names = [f.stem.replace("_pred_data", "") for f in files]
    res_dict = {
        problem_name: file_path for problem_name, file_path in zip(problem_names, files)
    }
    return res_dict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run prediction analysis for MOO-MAP experiments."
    )
    parser.add_argument(
        "--scenario_ids",
        type=int,
        nargs="*",
        help="IDs of the scenario to run (0-31). If not provided, all scenarios will be run.",
        default=None,
        required=False,
    )
    parser.add_argument(
        "--data_dir",
        type=str,
        help="Folder containing the prediction data files.",
        default=str(Path(__file__).parent / "results"),
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Folder to save the feature importance results.",
        default=str(Path(__file__).parent / "pred_results"),
    )
    parser.add_argument(
        "--file_ids",
        type=int,
        nargs="*",
        help="IDs of the data file to run (usually 0-22). If not provided, all files will be run.",
        default=None,
        required=False,
    )
    parser.add_argument(
        "--batch_id",
        type=int,
        nargs="*",
        help="IDs of the batch to run (combination of files and scenarios). If not provided, all combinations will be run.",
        default=None,
        required=False,
    )
    parser.add_argument(
        "--max_iter",
        type=int,
        default=200,
        help="Maximum number of iterations for the model.",
    )
    parser.add_argument(
        "--learning_rate", type=float, default=0.05, help="Learning rate for the model."
    )
    parser.add_argument(
        "--max_leaf_nodes",
        type=int,
        default=15,
        help="Maximum number of leaf nodes for the model.",
    )
    parser.add_argument(
        "--l2_regularization",
        type=float,
        default=1.0,
        help="L2 regularization strength for the model.",
    )
    parser.add_argument(
        "--n_splits", type=int, default=5, help="Number of splits for cross-validation."
    )
    parser.add_argument(
        "--n_repeats",
        type=int,
        default=3,
        help="Number of repeats for cross-validation.",
    )
    parser.add_argument(
        "--n_permutation_repeats",
        type=int,
        default=10,
        help="Number of repeats for permutation importance.",
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducibility."
    )
    args = parser.parse_args()
    return args


def get_config_dicts(args):
    scenario_list = generate_scenario_list()
    file_list = find_data_files(args.data_dir)
    scenario_ids = (
        args.scenario_ids
        if args.scenario_ids is not None
        else list(range(len(scenario_list)))
    )
    assert all(
        scenario_id < len(scenario_list) for scenario_id in scenario_ids
    ), f"Scenario IDs must be between 0 and {len(scenario_list)-1}"
    file_ids = (
        args.file_ids if args.file_ids is not None else list(range(len(file_list)))
    )
    assert all(
        file_id < len(file_list) for file_id in file_ids
    ), f"File IDs must be between 0 and {len(file_list)-1}"
    seed = (
        args.seed
        if args.seed is not None
        else int(np.random.SeedSequence().generate_state(1)[0])
    )

    batch_list = product(scenario_ids, file_ids)
    batch_ids = (
        args.batch_id
        if args.batch_id is not None
        else list(range(len(list(batch_list))))
    )

    config_dicts = []
    for batch_id in batch_ids:
        scenario_id, file_id = list(product(scenario_ids, file_ids))[batch_id]
        scenario = scenario_list[scenario_id]
        problem_name, file_path = list(file_list.items())[file_id]
        config_dict = {
            "file_path": file_path,
            "input_cols": scenario["input_cols"],
            "target_col": scenario["target_col"],
            "target_type": scenario["target_type"],
            "max_iter": args.max_iter,
            "learning_rate": args.learning_rate,
            "max_leaf_nodes": args.max_leaf_nodes,
            "l2_regularization": args.l2_regularization,
            "n_splits": args.n_splits,
            "n_repeats": args.n_repeats,
            "n_permutation_repeats": args.n_permutation_repeats,
            "random_state": seed,
            "problem_name": problem_name,
            "file_id": file_id,
            "scenario_id": scenario_id,
            "batch_id": batch_id,
        }
        config_dicts.append(config_dict)
    return config_dicts


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
