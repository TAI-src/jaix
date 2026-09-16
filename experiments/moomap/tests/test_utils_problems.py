import numpy as np
import pytest
from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)

from utils_cobi_configs import get_config
from utils_cobi_configs import names as cobi_names
from utils_problems import (
    cobi_problem_list,
    generate_problem_list,
    get_problem_info,
    get_problem_names,
    re_problem_list,
    get_cocoviz_problem_description,
)


def test_re_problem_list():
    problems = re_problem_list()
    assert isinstance(problems, list)
    assert len(problems) == 16  # There are 16 RE problems
    assert all(isinstance(p, REProblem) for p in problems)


@pytest.mark.parametrize("problem_ids", [[7, 8, 9], [10, 11], [12], [0]])
def test_re_problem_list_with_ids(problem_ids):
    if any(i < 7 or i >= 23 for i in problem_ids):
        with pytest.raises(AssertionError):
            re_problem_list(problem_ids)
        return
    problems = re_problem_list(problem_ids)
    assert isinstance(problems, list)
    assert len(problems) == len(problem_ids)
    assert all(isinstance(p, REProblem) for p in problems)

    problem_names = list(REProblem.problem_map.keys())
    for i, problem in zip(problem_ids, problems):
        expected_name = problem_names[i - 7]
        assert problem.problem_name == expected_name


def test_cobi_problem_list():
    problems = cobi_problem_list()
    assert isinstance(problems, list)
    assert len(problems) == 7  # There are 7 Cobi problems
    assert all(isinstance(p, CobiProblem) for p in problems)


@pytest.mark.parametrize("problem_ids", [[0, 1, 2], [3, 4], [5], [8]])
def test_cobi_problem_list_with_ids(problem_ids):
    if any(i < 0 or i >= 7 for i in problem_ids):
        with pytest.raises(AssertionError):
            cobi_problem_list(problem_ids)
        return
    problems = cobi_problem_list(problem_ids)
    assert isinstance(problems, list)
    assert len(problems) == len(problem_ids)
    assert all(isinstance(p, CobiProblem) for p in problems)

    for i, problem in zip(problem_ids, problems):
        expected_name = cobi_names[i]
        assert problem.name == expected_name


def test_generate_problem_list():
    problems = generate_problem_list()
    assert isinstance(problems, list)
    assert len(problems) == 23  # 7 Cobi + 16 RE problems
    assert all(isinstance(p, (CobiProblem, REProblem)) for p in problems)


@pytest.mark.parametrize("problem_ids", [[0, 1, 2], [7, 8], [0, 7, 15], [22], [23]])
def test_generate_problem_list_with_ids(problem_ids):
    if any(i < 0 or i >= 23 for i in problem_ids):
        with pytest.raises(AssertionError):
            generate_problem_list(problem_ids)
        return
    problems = generate_problem_list(problem_ids)
    assert isinstance(problems, list)
    assert len(problems) == len(problem_ids)
    assert all(isinstance(p, (CobiProblem, REProblem)) for p in problems)

    for i, problem in zip(problem_ids, problems):
        if 0 <= i < 7:
            expected_name = cobi_names[i]
            assert isinstance(problem, CobiProblem)
            assert problem.name == expected_name
        elif 7 <= i < 23:
            expected_name = list(REProblem.problem_map.keys())[i - 7]
            assert isinstance(problem, REProblem)
            assert problem.problem_name == expected_name


@pytest.mark.parametrize(
    "problem",
    [0, 7, REProblem(REProblemConfig(), 0), CobiProblem(get_config(0), inst=0)],
)
def test_get_problem_info(problem):
    info = get_problem_info(problem)
    assert isinstance(info, dict)
    assert "problem" in info
    assert "ideal_point" in info
    assert "nadir_point" in info
    assert "num_objectives" in info
    assert "num_variables" in info
    assert "lower_bounds" in info
    assert "upper_bounds" in info
    assert isinstance(info["problem"], str)
    assert isinstance(info["ideal_point"], np.ndarray)


def test_get_problem_names():
    names_dict = get_problem_names()
    assert isinstance(names_dict, dict)
    assert len(names_dict) == 23  # 7 Cobi + 16 RE problems
    for i in range(23):
        assert i in names_dict
        if 0 <= i < 7:
            assert names_dict[i] == cobi_names[i]
        elif 7 <= i < 23:
            expected_name = list(REProblem.problem_map.keys())[i - 7]
            assert names_dict[i] == expected_name


def test_get_cocoviz_problem_description():
    problem_names = get_problem_names()

    problem = 0  # CobiProblem
    desc = get_cocoviz_problem_description(problem)
    assert desc.name == problem_names[problem]
    assert desc.instance == 0
    assert desc.number_of_variables == 2
    assert desc.number_of_objectives == 2
    problem = 7  # REProblem
    desc = get_cocoviz_problem_description(problem)
    assert desc.name == f"REProblem_{problem_names[problem]}"
    assert desc.instance == 0
    assert desc.number_of_variables == 4
    assert desc.number_of_objectives == 2

    problem = REProblem(REProblemConfig(), 0)
    desc = get_cocoviz_problem_description(problem)
    assert desc.name == f"REProblem_{problem_names[7]}"
    assert desc.instance == 0
    assert desc.number_of_variables == 4
    assert desc.number_of_objectives == 2
