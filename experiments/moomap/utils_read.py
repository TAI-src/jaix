import json
import os
from collections import defaultdict
from pathlib import Path
import pandas as pd

from utils_problems import get_problem_info, get_problem_names
from plots_parallel_coordinate_plot import plot_parallel_coordinate_plot


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
            problem_idx = int(f.stem.split("_")[0])
            res_dict[problem_idx].append(f)
    return dict(res_dict)


def get_feature_importance_per_scenario(
    results_dir: str | Path,
    problem_ids: list[int] | None = None,
):
    data_files = find_data_files(
        results_dir, file_type_pattern="*_feat_imp.csv", problem_ids=problem_ids
    )
    data_file_by_scenario = defaultdict(list)
    for problem_id, files in data_files.items():
        for f in files:
            scenario_id = f.stem.split("_")[-3]
            assert scenario_id.startswith(
                "s"
            ), f"Scenario ID {scenario_id} does not start with 's'"
            data_file_by_scenario[scenario_id].append((problem_id, f))

    for scenario_id, file_tuple in data_file_by_scenario.items():
        scenario_df = pd.DataFrame()
        for problem_id, f in file_tuple:
            df = pd.read_csv(f)
            df["problem_id"] = problem_id
            scenario_df = pd.concat([scenario_df, df], ignore_index=True)
        scenario_df = scenario_df.rename(columns={"Unnamed: 0": "feature"})
        # average the feature importance over all problems for this scenario
        avg_vals = scenario_df.groupby("feature").mean().reset_index()
        avg_vals["problem_id"] = (
            100  # indicate that this is the average over all problems
        )
        scenario_df = pd.concat([scenario_df, avg_vals], ignore_index=True)
        for problem_id in list(scenario_df["problem_id"].unique()):
            problem_df = scenario_df[scenario_df["problem_id"] == problem_id]
            # plot feature importance as a parralel coordinates plot

            print(problem_df.columns)
            plot_parallel_coordinate_plot(
                problem_df,
                class_column="feature",
                features=[
                    "mutual_info",
                    "perm_importance",
                    "loo_cv_score",
                    "loo_score_drop",
                ],
                save_path=f"{results_dir}/feature_importance_{scenario_id}_problem_{problem_id}.pdf",
            )
