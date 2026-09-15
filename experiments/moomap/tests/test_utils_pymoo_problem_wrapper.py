from utils_pymoo_problem_wrapper import PymooProblemWrapper
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from pymoo.core.problem import ElementwiseProblem

from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.optimize import minimize
from pymoo.util.ref_dirs import get_reference_directions
from jaix.env.utils.archive.mo_archive import KeepDominated


def test_create_eval_archive():
    problem = REProblem(REProblemConfig(), inst=0)
    pymoo_problem = PymooProblemWrapper(problem)

    # Check if the archive is created correctly
    assert pymoo_problem.archive is not None
    assert pymoo_problem.archive.config.max_size is None
    assert pymoo_problem.archive.config.keep_dominated == KeepDominated.NONE
    assert pymoo_problem.archive.config.only_new_entries is False
    assert (
        pymoo_problem.archive.config.secondary_criterion_class.__name__
        == "ReferenceVectorDistanceScorer"
    )


def test_init_archive():
    problem = REProblem(REProblemConfig(), inst=0)
    pymoo_problem = PymooProblemWrapper(problem)

    # Check if the archive is initialized correctly
    assert pymoo_problem.archive is not None
    assert isinstance(pymoo_problem, ElementwiseProblem)


def test_with_pymoo():
    problem = REProblem(REProblemConfig(), inst=0)
    pymoo_problem = PymooProblemWrapper(problem)
    # create the reference directions to be used for the optimization
    ref_dirs = get_reference_directions(
        "das-dennis", problem.num_objectives, n_partitions=12
    )

    # create the algorithm object
    algorithm = NSGA3(pop_size=92, ref_dirs=ref_dirs)

    # execute the optimization
    res = minimize(pymoo_problem, algorithm, seed=1, termination=("n_gen", 2))
    assert len(res.opt) > 0
    assert pymoo_problem.archive.size > 0
