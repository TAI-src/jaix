import json
import os
from collections import defaultdict
from pathlib import Path

from utils_problems import generate_problem_list, get_problem_info, get_problem_names


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
    problem_list = generate_problem_list(problem_ids)
    problem_ids = (
        [i for i in range(len(problem_list))] if problem_ids is None else problem_ids
    )
    problem_infos = {
        i: get_problem_info(problem) for i, problem in zip(problem_ids, problem_list)
    }

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
    return config[config_type] if config_type in config else config


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
