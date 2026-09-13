from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)

from utils_cobi_configs import get_config
from utils_cobi_configs import names as cobi_names


def re_problem_list(problem_ids: list[int] | None = None):
    # RE problems are indexed from 7 to 22 (inclusive) in the combined problem list.
    assert problem_ids is None or all(
        isinstance(i, int) and 7 <= i < 23 for i in problem_ids
    ), "Problem IDs must be integers between 7 and 22"
    re_instances = list(range(7, 23)) if problem_ids is None else problem_ids
    # All non-constrained RE problems are included in the list. The constrained ones are excluded for now.
    problems = [REProblem(REProblemConfig(), i - 7) for i in re_instances]
    return problems


def cobi_problem_list(problem_ids: list[int] | None = None):
    assert problem_ids is None or all(
        isinstance(i, int) and 0 <= i < 7 for i in problem_ids
    ), "Problem IDs must be integers between 0 and 6"
    cobi_fun_ids = list(range(7)) if problem_ids is None else problem_ids
    cobi_configs = [get_config(func_id) for func_id in cobi_fun_ids]
    problems = [CobiProblem(config, inst=0) for config in cobi_configs]
    for problem_id in cobi_fun_ids:
        prob_idx = cobi_fun_ids.index(problem_id)
        problems[prob_idx].name = cobi_names[problem_id]

    return problems


def generate_problem_list(problem_ids: list[int] | None = None):
    assert problem_ids is None or all(
        isinstance(i, int) and 0 <= i < 23 for i in problem_ids
    ), "Problem IDs must be integers between 0 and 22"
    if problem_ids is not None:
        cobi_ids = [i for i in problem_ids if 0 <= i < 7]
        re_ids = [i for i in problem_ids if 7 <= i < 23]
    else:
        cobi_ids = None
        re_ids = None
    cobi_problems = cobi_problem_list(cobi_ids)
    re_problems = re_problem_list(re_ids)
    return cobi_problems + re_problems


def get_problem_info(problem: CobiProblem | REProblem | int) -> dict:
    if isinstance(problem, int):
        problem_list = generate_problem_list(problem_ids=[problem])
        problem = problem_list[0]
    return {
        "problem": str(problem),
        "ideal_point": problem.ideal_point,
        "nadir_point": problem.nadir_point,
        "num_objectives": problem.num_objectives,
        "num_variables": problem.dimension,
        "lower_bounds": problem.lower_bounds,
        "upper_bounds": problem.upper_bounds,
    }


def get_problem_names(problem_ids: list[int] | None = None) -> dict[int, str]:
    cobi_names_dict = {i: name for i, name in enumerate(cobi_names)}
    re_names_dict = {i + 7: name for i, name in enumerate(REProblem.problem_map.keys())}
    full_names_dict = {**cobi_names_dict, **re_names_dict}
    if problem_ids is None:
        problem_ids = list(range(23))
    assert all(
        isinstance(i, int) and 0 <= i < 23 for i in problem_ids
    ), "Problem IDs must be integers between 0 and 22"
    return {i: full_names_dict[i] for i in problem_ids if i in full_names_dict}
