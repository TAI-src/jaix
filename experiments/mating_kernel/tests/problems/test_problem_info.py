import pytest
from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)

from mating_kernel.problems.cobi_configs import get_config
from mating_kernel.problems.problem_info import ProblemInfo

cobi_problem = CobiProblem(get_config(0), inst=0)
re_problem = REProblem(REProblemConfig(), 0)


@pytest.mark.parametrize("problem", [cobi_problem, re_problem])
def test_problem_info(problem):
    problem_info = ProblemInfo(problem)
    assert problem_info.problem == str(problem)
    assert hasattr(problem_info, "ideal_point")

    assert problem_info.num_objectives == problem.num_objectives
    assert problem_info.num_variables == problem.dimension
    assert hasattr(problem_info, "uuid")
