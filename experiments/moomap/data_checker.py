import csv
import json
import os
from collections import defaultdict
from pathlib import Path

from utils_read import get_nsga3x_results

results_dir = f"{os.path.dirname(os.path.abspath(__file__))}/results"
EXPECTED_NUM_PROBLEMS = 23
EXPECTED_NUM_RUNS = 30

results = get_nsga3x_results(results_dir=results_dir)

# Check that all problem ids are present
assert set(results.keys()) == set(
    range(EXPECTED_NUM_PROBLEMS)
), f"Missing problem ids {set(range(EXPECTED_NUM_PROBLEMS)) - set(results.keys())}"
for problem_id, problem_dict in results.items():
    print(
        f"Checking results for problem: {problem_dict["problem"]} with index {problem_id}/{len(results)}"
    )
    runs = [run_key for run_key in problem_dict.keys() if run_key.startswith("r_")]
    if len(runs) != EXPECTED_NUM_RUNS:
        print(
            f"Unexpected number of runs for problem: {problem_dict['problem']}. Expected {EXPECTED_NUM_RUNS}, got {len(runs)}"
        )
    seeds = []
    for rid in runs:
        run_dict = problem_dict[rid]
        seeds.append(run_dict["config"]["seed"])
        n_gen = run_dict["config"]["num_generations"]
        n_prefill = run_dict["config"]["num_prefill_samples"]
        n_offspring = run_dict["config"]["num_offspring"]
        assert (
            n_prefill == n_offspring
        ), f"Unexpected prefill samples for problem: {problem_dict['problem']}. Expected {n_offspring}, got {n_prefill}"
        expected_num_results = n_gen * n_offspring + 1  # plus 1 for the header
        # Count the number of lines in the result file
        res_file = run_dict["result_file"]
        with open(res_file, "r") as rf:
            reader = csv.reader(rf)
            num_lines = sum(1 for line in reader)
            if num_lines != expected_num_results:
                print(
                    f"Unexpected number of results in file: {res_file}. Expected {expected_num_results}, got {num_lines}"
                )

    # Check for duplicate seeds
    if len(seeds) != len(set(seeds)):
        print(
            f"Duplicate seeds found for problem: {problem_dict['problem']}. Seeds: {seeds}"
        )
