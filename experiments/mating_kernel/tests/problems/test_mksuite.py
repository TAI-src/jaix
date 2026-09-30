import pytest
from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
)

from mating_kernel.problems.mk_suite import MKSuite, MKSuiteConfig
from mating_kernel.problems.problem_info import ProblemInfo


def test_re_problem_list():
    problems = MKSuite.re_problem_list(constrained=False)
    assert isinstance(problems, list)
    assert len(problems) == 16
    assert all(isinstance(p, REProblem) for p in problems)


def test_cobi_problem_list():
    problems = MKSuite.cobi_problem_list()
    assert isinstance(problems, list)
    assert len(problems) == 7
    assert all(hasattr(p, "name") for p in problems)
    assert problems[0].name == "cobi_lin"
    # Quick integration test to ensure that the name is override is passed on
    pinfo = ProblemInfo(problems[0])
    assert pinfo.uuid == "cobi_lin"
    assert all(isinstance(p, CobiProblem) for p in problems)


@pytest.mark.parametrize(
    "cobi, re, constrained, num_objectives",
    [
        (True, True, False, None),
        (True, False, False, None),
        (False, True, False, None),
        (True, True, True, None),
        (True, True, False, [2]),
        (True, True, False, [2, 3]),
    ],
)
def test_generate_problem_list(cobi, re, constrained, num_objectives):
    problems = MKSuite.generate_problem_list(cobi, re, constrained, num_objectives)
    assert isinstance(problems, list)
    if cobi:
        if not re:
            assert all(isinstance(p, CobiProblem) for p in problems)
        else:
            assert any(isinstance(p, CobiProblem) for p in problems)
    if re:
        if not cobi:
            assert all(isinstance(p, REProblem) for p in problems)
        else:
            assert any(isinstance(p, REProblem) for p in problems)
    if num_objectives is not None:
        assert all(p.num_objectives in num_objectives for p in problems)


def test_init():
    config = MKSuiteConfig(cobi=True, re=True, constrained=False, num_objectives=None)
    mk_suite = MKSuite(config)
    assert hasattr(mk_suite, "problems")
    assert hasattr(mk_suite, "problem_info_list")
    assert hasattr(mk_suite, "problem_id_map")
    assert len(mk_suite.problem_id_map) == len(mk_suite.problems)
    assert len(mk_suite.problems) == 23  # 7 Cobi + 16 RE problems
    assert mk_suite.problem_id_map["cobi_lin"]["info"].uuid == "cobi_lin"
    assert isinstance(mk_suite.problem_id_map["cobi_lin"]["problem"], CobiProblem)
