from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.optimize import minimize
from pymoo.problems import get_problem
from pymoo.util.ref_dirs import get_reference_directions

from utils_archive_stats_callback import ArchiveStatsCallback
from utils_pymoo_problem_wrapper import PymooProblemWrapper


def test_archive_stats_callback():
    # create the reference directions to be used for the optimization
    ref_dirs = get_reference_directions("das-dennis", 3, n_partitions=12)

    # create the algorithm object
    algorithm = NSGA3(pop_size=92, ref_dirs=ref_dirs)

    archive = PymooProblemWrapper._create_eval_archive(
        REProblem(REProblemConfig(), inst=0)
    )
    callback = ArchiveStatsCallback(archive=archive)

    # execute the optimization
    minimize(
        get_problem("dtlz1"),
        algorithm,
        seed=1,
        termination=("n_gen", 5),
        callback=callback,
    )
    assert len(callback.data["archive_stats"]) == 5
