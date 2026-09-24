import json
import os
from collections import defaultdict
from pathlib import Path
import pandas as pd

from utils_problems import get_problem_info, get_problem_names
from plots_parallel_coordinate_plot import plot_pcp


def get_nsga3x_results(
    results_dir: str | Path | None = None, problem_ids: list[int] | None = None
) -> dict:
    """
    Get the results of the NSGA3x experiments.
    """
    if results_dir is None:
        results_dir = f"{os.path.dirname(os.path.abspath(__file__))}/results"
    result_files = find_data_files(
        results_dir, file_type_pattern="results_*.csv", problem_ids=problem_ids
    )
    config_files = find_data_files(
        results_dir, file_type_pattern="config_*.json", problem_ids=problem_ids
    )
    if problem_ids is None:
        problem_ids = list(get_problem_names().keys())
    problem_infos = {i: get_problem_info(i) for i in problem_ids}

    for problem_id, problem_info in problem_infos.items():
        resf = result_files.get(problem_id, [])
        conf = config_files.get(problem_id, [])
        assert len(resf) == len(
            conf
        ), f"Number of result files and config files do not match for problem {problem_id}"
        for res_file, config_file in zip(resf, conf):
            run_id = res_file.parent.name
            assert (
                run_id == config_file.parent.name
            ), f"Run ID mismatch: {run_id} != {config_file.parent.name}"
            config = get_config_dict(config_file)
            problem_info[run_id] = {
                "result_file": res_file,
                "config_file": config_file,
                "config": config,
            }
    return problem_infos


def get_pred_overview_results(
    results_dir: str | Path, problem_ids: list[int] | None = None
) -> pd.DataFrame:
    config_files = find_data_files(
        results_dir, file_type_pattern="*_config.json", problem_ids=problem_ids
    )
    data = []
    for _, files in config_files.items():
        for config_file in files:
            cfg = get_config_dict(config_file)
            data.append(
                {
                    "problem_id": cfg.get("problem_id", cfg.get("problem_name")),
                    "scenario_id": cfg.get("scenario_id"),
                    "cv_score_mean": cfg.get("cv_score_mean"),
                    "batch_id": cfg.get("batch_id"),
                }
            )

    df = pd.DataFrame(data)
    return df


def get_config_dict(config_file, config_type="NSGA3ExperimentConfig") -> dict:
    """
    Get the config dict from a config file.
    """
    with open(config_file, "r") as f:
        config = json.load(f)
    return config.get(config_type, config)


def find_data_files(
    folder: str | Path,
    file_type_pattern: str = "",
    problem_ids: list[int] | None = None,
) -> dict[int, list[Path]]:
    path = Path(folder)
    files = list(path.rglob(file_type_pattern))

    res_dict = defaultdict(list)
    for f in files:
        name_not_found = True
        for problem_idx, name in get_problem_names(problem_ids).items():
            if name in f.name:
                res_dict[problem_idx].append(f)
                name_not_found = False
                break
        if name_not_found:
            try:
                problem_idx = int(f.stem.split("_")[0])
                res_dict[problem_idx].append(f)
            except ValueError:
                print(
                    f"Could not find problem name for file {f.name}. Please check the file name and ensure it contains a valid problem name or index."
                )
    return dict(res_dict)


def get_feature_importance_per_scenario(
    results_dir: str | Path,
    problem_ids: list[int] | None = None,
):
    data_files = find_data_files(
        results_dir, file_type_pattern="*_feat_imp.csv", problem_ids=problem_ids
    )
    features = [
        "mutual_info",
        "perm_importance",
        "loo_score_drop",
    ]

    data_file_by_scenario = defaultdict(list)
    for problem_id, files in data_files.items():
        for f in files:
            scenario_id = f.stem.split("_")[-3]
            assert scenario_id.startswith(
                "s"
            ), f"Scenario ID {scenario_id} does not start with 's'"
            data_file_by_scenario[scenario_id].append((problem_id, f))

    average_dfs = []
    for scenario_id, file_tuple in data_file_by_scenario.items():
        scenario_df = pd.DataFrame()
        for problem_id, f in file_tuple:
            df = pd.read_csv(f)
            df["problem_id"] = problem_id
            # only keep features of interest
            df = df[["Unnamed: 0"] + features + ["problem_id"]]
            scenario_df = pd.concat([scenario_df, df], ignore_index=True)
        scenario_df = scenario_df.rename(columns={"Unnamed: 0": "feature"})
        # check if all the loo_score_drop values are negative. If so, this is old data that swapped the values around, so we need to swap them back
        if (scenario_df["loo_score_drop"] < 0).all():
            scenario_df["loo_score_drop"] = -scenario_df["loo_score_drop"]
        # average the feature importance over all problems for this scenario
        avg_vals = scenario_df.groupby("feature").mean().reset_index()
        avg_vals["problem_id"] = "avg"
        scenario_df = pd.concat([scenario_df, avg_vals], ignore_index=True)
        avg_vals["scenario_id"] = scenario_id
        average_dfs.append(avg_vals)

        for problem_id in list(scenario_df["problem_id"].unique()):

            problem_df = scenario_df[
                (scenario_df["problem_id"] == problem_id)
                | (scenario_df["problem_id"] == "avg")
            ]

            plot_df = (
                problem_df.set_index(["problem_id", "feature"])
                .T.stack(level=0)
                .reset_index()
                .rename(columns={"level_0": "metric", "level_1": "problem_id"})
            )

            # plot feature importance as a parralel coordinates plot
            plot_pcp(
                plot_df,
                line_col_name="metric",
                class_col_name="problem_id",
                linestyles={"avg": "--"},
                save_path=f"{results_dir}/feat_{scenario_id}_p{problem_id}.pdf",
            )
    average_df = pd.concat(average_dfs, ignore_index=True)
    # drop scenario_id and problem_id columns
    average_df = average_df.drop(columns=["scenario_id", "problem_id"])
    avg_all = (
        average_df.groupby("feature")
        .mean()
        .reset_index()
        .set_index("feature")
        .T.reset_index()
        .rename(columns={"index": "metric"})
    )
    plot_pcp(
        avg_all,
        line_col_name="metric",
        save_path=f"{results_dir}/feat_avg_all.pdf",
    )
