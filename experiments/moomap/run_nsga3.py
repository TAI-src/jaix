import os

import pandas as pd
from jaix.env.utils.mo_sizing import get_ref_dirs
from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import REProblem
from pymoo.algorithms.moo.nsga3 import NSGA3, ReferenceDirectionSurvival
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize

from config_run_nsga3 import get_batches, parse_args
from utils_archive_stats_callback import ArchiveStatsCallback
from utils_nsga3_norm import StaticReferenceDirectionSurvival
from utils_pymoo_problem_wrapper import PymooProblemWrapper


def run_algorithm(
    n_gen: int,
    seed: int,
    out_dir: str,
    problem: REProblem | CobiProblem,
    static_ref: bool = True,
    verbose: bool = False,
    algorithm_name: str = "nsga3",
) -> str:
    ref_dirs = get_ref_dirs(problem.num_objectives, "original")

    # create the algorithm object
    if static_ref:
        survival = StaticReferenceDirectionSurvival(
            ref_dirs, problem.ideal_point, problem.nadir_point
        )
    else:
        survival = ReferenceDirectionSurvival(ref_dirs)
    if algorithm_name == "nsga2":
        algorithm = NSGA2(pop_size=len(ref_dirs))
        static_ref = False  # NSGA2 does not use reference directions
    elif algorithm_name == "nsga3":
        algorithm = NSGA3(
            pop_size=len(ref_dirs),
            ref_dirs=ref_dirs,
            survival=survival,
        )
    else:
        raise ValueError(
            f"Invalid algorithm name: {algorithm_name}. Choose 'nsga2' or 'nsga3'."
        )

    pymoo_problem = PymooProblemWrapper(problem)
    callback = ArchiveStatsCallback(archive=pymoo_problem.archive)

    # execute the optimization
    minimize(
        pymoo_problem,
        algorithm,
        seed=seed,
        termination=("n_gen", n_gen),
        callback=callback,
        verbose=verbose,
    )
    file_name = (
        f"{algorithm_name}_{problem!s}_s{seed}{"_fixed" if static_ref else ""}.csv"
    )
    file_path = f"{out_dir}/{file_name}"
    df = pd.DataFrame(callback.data["archive_stats"])
    df.to_csv(file_path, index=False)
    return file_path


def run(args):
    os.makedirs(args.out_dir, exist_ok=True)
    batch_list, num_batches = get_batches(
        batch_id=args.batch_ids,
        n_runs=args.n_runs,
        problem_ids=args.problem_ids,
        static_ref=args.static_ref,
        seed=args.seed,
    )
    files = []
    for batch in batch_list:
        print("Running batch:", batch["id"], "out of", num_batches)
        out_file = run_algorithm(
            n_gen=args.n_gen,
            seed=batch["seed"],
            out_dir=args.out_dir,
            static_ref=batch["static_ref"],
            problem=batch["problem"],
            verbose=args.verbose,
            algorithm_name=args.algorithm,
        )
        files.append(out_file)
    return files


if __name__ == "__main__":
    args = parse_args()
    run(args)
