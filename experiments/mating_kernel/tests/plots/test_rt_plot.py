import pytest

from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)

from mating_kernel.problems.cobi_configs import get_config
from mating_kernel.problems.problem_info import ProblemInfo
from mating_kernel.plots.rt_plot import get_description

cobi_problem = CobiProblem(get_config(0), inst=0)
re_problem = REProblem(REProblemConfig(), 0)


test_data_regular = [
    {
        "a_coverage": 0.5,
        "a_score": 0.6,
        "a_fevals": 200,
        "seed": 1,
    },
    {
        "a_coverage": 0.5,
        "a_score": 0.6,
        "a_fevals": 200,
        "seed": 1,
    },
    {
        "a_coverage": 0.6,
        "a_score": 0.7,
        "a_fevals": 300,
        "seed": 1,
    },
    {
        "a_coverage": 0.6,
        "a_score": 0.7,
        "a_fevals": 300,
        "seed": 1,
    },
    {
        "a_coverage": 0.7,
        "a_score": 0.8,
        "a_fevals": 200,
        "seed": 2,
    },
]


@pytest.mark.parametrize("problem", [cobi_problem, re_problem])
def test_get_description(problem):
    problem_info = ProblemInfo(problem)
    get_description(problem_info)  # Ensure it runs without errors
