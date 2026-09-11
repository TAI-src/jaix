import os
from collections import defaultdict
from pathlib import Path
import json

from cobi_config_generator import old_names
from nsga3_experiment import NSGA3ExperimentConfig


def build_results_dict(results_dir: str | None = None) -> dict:
    """
    Build a dictionary of results from the logs directory.
    """
    if results_dir is None:
        results_dir = f"{os.path.dirname(os.path.abspath(__file__))}/results"
    res_dict = defaultdict(dict)
    # Get all csv files in the logs directory
    csv_files = list(Path(results_dir).rglob("*.csv"))
    for file in csv_files:
        problem_str = file.name.replace("results_", "").replace(".csv", "")
        tmp_dict = {"file": file}

        # get example config to know number of niches (dirty)
        # remove results_ prefix and .csv suffix from problem name to get config name
        config_file = file.parent / f"config_{problem_str}.json"
        # load the config to get the number of niches

        with open(config_file, "r") as f:
            config = json.load(f)
            tmp_dict["config"] = config["NSGA3ExperimentConfig"]

        res_dict[problem_str][file.parent.name] = tmp_dict

    res_dict = dict(res_dict)  # convert to normal dict
    expected_problem_list = NSGA3ExperimentConfig.generate_problem_list()
    problem_strs = [str(problem) for problem in expected_problem_list]
    for problem_str, problem in zip(problem_strs, expected_problem_list):
        if problem_str.startswith("cobi_"):
            # translate old names to new names
            problem_str = old_names.get(problem_str, problem_str)
        if problem_str not in res_dict:
            print(f"Missing results for problem: {problem_str}")
        else:
            # add info about the problem to the dict
            res_dict[problem_str]["ideal_point"] = problem.ideal_point
            res_dict[problem_str]["nadir_point"] = problem.nadir_point
            res_dict[problem_str]["num_objectives"] = problem.num_objectives
            res_dict[problem_str]["num_variables"] = problem.dimension
            res_dict[problem_str]["lower_bounds"] = problem.lower_bounds
            res_dict[problem_str]["upper_bounds"] = problem.upper_bounds

    return res_dict
