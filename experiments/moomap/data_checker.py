import csv
import os
import argparse
from utils_read import get_nsga3x_results, get_pred_overview_results
from utils_problems import get_problem_names

EXPECTED_NUM_PROBLEMS = 23
EXPECTED_NUM_RUNS = 30
EXPECTED_NUM_SCENARIOS = 32


def check_experiment_results(
    results_dir: str,
    expected_num_problems: int = EXPECTED_NUM_PROBLEMS,
    expected_num_runs: int = EXPECTED_NUM_RUNS,
):
    results = get_nsga3x_results(results_dir=results_dir)

    # Check that all problem ids are present
    assert set(results.keys()) == set(
        range(expected_num_problems)
    ), f"Missing problem ids {set(range(expected_num_problems)) - set(results.keys())}"
    for problem_id, problem_dict in results.items():
        print(
            f"Checking results for problem: {problem_dict["problem"]} with index {problem_id}/{len(results)}"
        )
        runs = [run_key for run_key in problem_dict if run_key.startswith("r_")]
        if len(runs) != expected_num_runs:
            print(
                f"Unexpected number of runs for problem: {problem_dict['problem']}. Expected {expected_num_runs}, got {len(runs)}"
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


def check_pred_results(
    results_dir: str,
    expected_num_problems: int = EXPECTED_NUM_PROBLEMS,
    expected_num_runs: int = EXPECTED_NUM_SCENARIOS,
):

    results = get_pred_overview_results(results_dir=results_dir)
    problem_ids = [int(pid) for pid in results["problem_id"].unique()]

    # Check that all problem ids are present
    assert set(problem_ids) == set(
        range(expected_num_problems)
    ), f"Missing problem ids {set(range(expected_num_problems)) - set(problem_ids)}"

    problem_names = get_problem_names(problem_ids=problem_ids)
    for problem_id, problem_name in problem_names.items():
        problem_results = results[results["problem_id"] == problem_id]
        # get the scenario ids and check for duplicates
        scenario_ids = problem_results["scenario_id"].tolist()
        if len(scenario_ids) != len(set(scenario_ids)):
            # print the duplicate scenario ids
            dupes = problem_results["scenario_id"].value_counts().loc[lambda x: x > 1]
            print(
                f"Duplicate scenario ids found for {problem_id}: {problem_name}. Scenario ids: {dupes.index.tolist()}"
            )

        # Check which scenarios are missing
        missing_scenarios = set(range(expected_num_runs)) - set(scenario_ids)
        if missing_scenarios:
            print(
                f"Missing scenario ids for {problem_id}: {problem_name}. Missing scenario ids: {missing_scenarios}"
            )
    batch_ids = results["batch_id"].unique()
    expected_batch_ids = list(range(expected_num_problems * expected_num_runs))
    missing_batches = set(expected_batch_ids) - set(batch_ids)
    if missing_batches:
        print(f"Missing batch ids: {missing_batches}")


def parse_args():

    parser = argparse.ArgumentParser(
        description="Check the results of the NSGA3X experiments."
    )
    parser.add_argument(
        "--results_dir",
        type=str,
        help="Directory containing the results of the experiments.",
    )
    parser.add_argument(
        "--expected_num_problems",
        type=int,
        default=EXPECTED_NUM_PROBLEMS,
        help="Expected number of problems.",
    )
    parser.add_argument(
        "--expected_num_runs",
        type=int,
        default=EXPECTED_NUM_RUNS,
        help="Expected number of runs per problem.",
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["experiment", "pred"],
        default="experiment",
        help="Mode of operation: 'experiment' to check NSGA3X experiment results, 'pred' to check prediction results.",
    )
    parser.add_argument(
        "--expected_num_scenarios",
        type=int,
        default=EXPECTED_NUM_SCENARIOS,
        help="Expected number of scenarios for prediction results.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.mode == "experiment":
        check_experiment_results(
            results_dir=args.results_dir,
            expected_num_problems=args.expected_num_problems,
            expected_num_runs=args.expected_num_runs,
        )
    elif args.mode == "pred":
        check_pred_results(
            results_dir=args.results_dir,
            expected_num_problems=args.expected_num_problems,
            expected_num_runs=args.expected_num_scenarios,
        )


if __name__ == "__main__":
    main()
