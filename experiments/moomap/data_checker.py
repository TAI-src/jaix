import csv
import json
import os
from collections import defaultdict
from pathlib import Path

from utils_problems import get_problem_names

logs_dir = f"{os.path.dirname(os.path.abspath(__file__))}/logs"
results_dir = f"{os.path.dirname(os.path.abspath(__file__))}/results"

# Get all csv files in the logs directory
csv_files = list(Path(results_dir).rglob("*.csv"))

res_dict = defaultdict(list)
for file in csv_files:
    res_dict[file.name].append(file)

problem_dict = get_problem_names()


problems = list(res_dict.keys())
for i, problem in problem_dict.items():
    print(f"Checking results for problem: {problem} with index {i}/{len(problem_dict)}")
    csv_name = f"results_{problem}.csv"
    if csv_name not in problems:
        print(f"Missing results for problem: {problem}")
    else:
        # Check the number of results are as expected
        config_name = f"config_{problem}.json"
        res_files = res_dict[csv_name]
        if len(res_files) != 30:
            print(
                f"Unexpected number of results for problem: {problem}. Expected 30, got {len(res_files)}"
            )

        seeds = []
        for res_file in res_files:
            # get the expected number of results from the corresponding config file
            config_file_path = res_file.parent / config_name
            if not config_file_path.exists():
                print(
                    f"Missing config file for problem: {problem} at {config_file_path}"
                )
                continue
            with open(config_file_path, "r") as f:
                config = json.load(f)
                n_gen = config["NSGA3ExperimentConfig"]["num_generations"]
                n_prefil = config["NSGA3ExperimentConfig"]["num_prefill_samples"]
                n_offspring = config["NSGA3ExperimentConfig"]["num_offspring"]
                seeds.append(config["NSGA3ExperimentConfig"]["seed"])
                expected_num_results = n_gen * n_offspring + 1  # plus 1 for the header
                # Count the number of lines in the result file
                with open(res_file, "r") as rf:
                    reader = csv.reader(rf)
                    num_lines = sum(1 for line in reader)
                    if num_lines != expected_num_results:
                        print(
                            f"Unexpected number of results in file: {res_file}. Expected {expected_num_results}, got {num_lines}"
                        )
        # Check if all seeds are unique
        if len(seeds) != len(set(seeds)):
            print(f"Duplicate seeds found for problem: {problem}. Seeds: {seeds}")
