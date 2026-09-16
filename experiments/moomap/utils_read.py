import json
import os
from collections import defaultdict
from pathlib import Path

from utils_problems import (
    get_problem_info,
    get_problem_names,
    get_cocoviz_problem_description,
)
from cocoviz import ProblemDescription, Result, ResultSet, Indicator, rtpplot
import pandas as pd


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


def get_config_dict(config_file, config_type="NSGA3ExperimentConfig") -> dict:
    """
    Get the config dict from a config file.
    """
    with open(config_file, "r") as f:
        config = json.load(f)
    return config.get(config_type, config)


def read_perf_results(
    results_dir: str | Path,
    algorithm_names: list[str],
    problem_ids: list[int] | None = None,
):

    all_results = defaultdict(dict)
    for algorithm in algorithm_names:
        res_dict = find_data_files(
            results_dir, file_type_pattern=f"{algorithm}_*.csv", problem_ids=problem_ids
        )
        for problem_id, files in res_dict.items():
            all_results[problem_id][algorithm] = files

    indicator_specs = [
        ("coverage", True),
        ("size", True),
        ("score", True),
        ("avg_dist_to_niches", False),
        ("avg_niche_count", True),
        ("filled_niches", True),
        ("avg_dist_to_ideal", False),
        ("niche_perf_avg", False),
    ]

    results = ResultSet()
    for problem_id, prob_files in all_results.items():
        problem_desc = get_cocoviz_problem_description(problem_id)
        for algo, files in prob_files.items():
            for file in files:
                df = pd.read_csv(file)
                result = Result(
                    algorithm=algo,
                    problem=problem_desc,
                    data=df,
                    fevals_column="fevals",
                )
                results.append(result)

    for problem, by_prob in results.by_problem_name():
        for ind_col, lib in indicator_specs:
            ind = Indicator(ind_col, larger_is_better=lib)
            ax = rtpplot(results, ind)
            # save the plot to a file
            file_name = f"{results_dir}/perf_{problem_id}_{ind_col}.pdf"
            ax.figure.savefig(file_name, dpi=300)


def find_data_files(
    folder: str | Path,
    file_type_pattern: str = "",
    problem_ids: list[int] | None = None,
) -> dict[int, list[Path]]:
    path = Path(folder)
    files = list(path.rglob(file_type_pattern))

    res_dict = defaultdict(list)
    for f in files:
        for problem_idx, name in get_problem_names(problem_ids).items():
            if name in f.name:
                res_dict[problem_idx].append(f)
                break
    return dict(res_dict)
