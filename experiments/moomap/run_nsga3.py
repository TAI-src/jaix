from pymoo.algorithms.moo.nsga3 import NSGA3, ReferenceDirectionSurvival
from utils_problems import generate_problem_list
from jaix.env.utils.mo_sizing import get_ref_dirs
from utils_pymoo_problem_wrapper import PymooProblemWrapper

from pymoo.optimize import minimize
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.visualization.scatter import Scatter

from utils_nsga3_norm import StaticReferenceDirectionSurvival
from utils_archive_stats_callback import ArchiveStatsCallback


def run_algorithm(
    n_gen: int, seed: int, static_ref: bool = True, problem_ids: list[int] | None = None
):

    for problem in generate_problem_list(problem_ids=problem_ids):
        ref_dirs = get_ref_dirs(problem.num_objectives, "original")

        # create the algorithm object
        if static_ref:
            survival = StaticReferenceDirectionSurvival(
                ref_dirs, problem.ideal_point, problem.nadir_point
            )
        else:
            survival = ReferenceDirectionSurvival(ref_dirs)
        algorithm = NSGA3(pop_size=len(ref_dirs), ref_dirs=ref_dirs, survival=survival)

        pymoo_problem = PymooProblemWrapper(problem)
        callback = ArchiveStatsCallback(archive=pymoo_problem.archive)

        # execute the optimization
        res = minimize(
            pymoo_problem,
            algorithm,
            seed=seed,
            termination=("n_gen", n_gen),
            callback=callback,
        )
        print(callback.data["archive_stats"])


if __name__ == "__main__":
    run_algorithm(n_gen=2, seed=42, static_ref=True, problem_ids=[0])
