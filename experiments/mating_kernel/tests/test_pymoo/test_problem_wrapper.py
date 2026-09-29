from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.core.problem import ElementwiseProblem
from pymoo.optimize import minimize
from pymoo.util.ref_dirs import get_reference_directions

from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper


def test_init():
    problem = REProblem(REProblemConfig(), inst=0)
    pymoo_problem = PymooProblemWrapper(problem)

    # Check if the archive is initialized correctly
    assert isinstance(pymoo_problem, ElementwiseProblem)
    assert len(pymoo_problem.records) == 0


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
    assert (
        len(pymoo_problem.records) == 2 * 92
    )  # 2 generations, 92 individuals per generation
    records = pymoo_problem.retrieve_records()
    assert len(records) == 2 * 92
    assert len(pymoo_problem.records) == 0


def test_without_recording():
    problem = REProblem(REProblemConfig(), inst=0)
    pymoo_problem = PymooProblemWrapper(problem, record=False)
    # create the reference directions to be used for the optimization
    ref_dirs = get_reference_directions(
        "das-dennis", problem.num_objectives, n_partitions=12
    )

    # create the algorithm object
    algorithm = NSGA3(pop_size=92, ref_dirs=ref_dirs)

    # execute the optimization
    res = minimize(pymoo_problem, algorithm, seed=1, termination=("n_gen", 2))
    assert len(res.opt) > 0
    assert len(pymoo_problem.records) == 0  # No records should be stored
