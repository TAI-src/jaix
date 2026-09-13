import argparse


def parse_args():

    parser = argparse.ArgumentParser(
        description="Postprocess results from NSGA3 experiments."
    )
    parser.add_argument(
        "--results_dir",
        type=str,
        default="results",
        help="Directory containing the results.",
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        default="postprocess_results",
        help="Directory to save the postprocessed results.",
    )
    parser.add_argument(
        "--skip_plots",
        action="store_true",
        help="Skip plotting the results.",
    )
    parser.add_argument(
        "--problem_ids",
        type=int,
        nargs="*",
        default=None,
        help="List of problem IDs to process. If not provided, all problems will be processed.",
    )
    return parser.parse_args()
