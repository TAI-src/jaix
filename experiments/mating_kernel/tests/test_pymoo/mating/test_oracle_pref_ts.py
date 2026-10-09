from mating_kernel.pymoo.mating.oracle_pref_ts import OraclePrefTournamentSelection
from .test_random_pref_ts import dummy_comp, create_pop
from pymoo.problems import get_problem
from copy import deepcopy
from unittest.mock import Mock

from mating_kernel.problems.mo_tracking import make_tracked
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper


def test_eval_mate():
    # Create a dummy problem and population
    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)
    pop_cpy = deepcopy(pop)  # Make a copy of the original population

    # Create an instance of OraclePrefTournamentSelection
    selection = OraclePrefTournamentSelection(
        func_comp=dummy_comp, num_oracle_offspring=5
    )

    # Select two parents from the population
    parent_idx = [0]  # Indices of the parents in the population
    mate_idx = 2  # Index of the mate in the population

    # Evaluate the mate
    survival_rate = selection.eval_mate(parent_idx, mate_idx, problem, pop)

    # Check that the survival rate is between 0 and 1
    assert 0 <= survival_rate <= 1
    assert survival_rate * 5 == round(
        survival_rate * 5
    )  # Because we have 5 offsprings, the survival rate should be a multiple of 1/5

    # Check that the original population has not been modified
    assert len(pop) == 10
    for ind_original, ind_copied in zip(pop, pop_cpy):
        assert ind_original.F.tolist() == ind_copied.F.tolist()
        assert ind_original.data == ind_copied.data


# Test with pymoo wrapper and archive to ensure the additional evals do not get added to the archive
def test_eval_mate_with_archive():
    tracked_REProblem = make_tracked(REProblem)
    static_problem = tracked_REProblem(REProblemConfig(), inst=0)
    problem = PymooProblemWrapper(static_problem)

    pop = create_pop(problem, size=10)
    pop_cpy = deepcopy(pop)  # Make a copy of the original population
    len_records = len(problem.records)  # Store the original length of records
    stats = (
        problem.static_problem.get_archive_stats()
    )  # Ensure the archive is initialized
    assert stats["size"] > 0  # Previous evaluations

    # Create an instance of OraclePrefTournamentSelection
    selection = OraclePrefTournamentSelection(
        func_comp=dummy_comp, num_oracle_offspring=5
    )

    # Select two parents from the population
    parent_idx = [0]  # Indices of the parents in the population
    mate_idx = 2  # Index of the mate in the population

    # Evaluate the mate
    survival_rate = selection.eval_mate(parent_idx, mate_idx, problem, pop)

    # Ensure that the archive and records have not been modified
    assert len(problem.records) == len_records
    assert problem.static_problem.get_archive_stats() == stats

    # Check that the survival rate is between 0 and 1
    assert 0 <= survival_rate <= 1

    # Check that the original population has not been modified
    assert len(pop) == 10
    for ind_original, ind_copied in zip(pop, pop_cpy):
        assert ind_original.F.tolist() == ind_copied.F.tolist()
        assert ind_original.data == ind_copied.data


def test_select_mate():
    # Create a dummy problem and population
    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)

    # Create an instance of OraclePrefTournamentSelection
    selection = OraclePrefTournamentSelection(
        func_comp=dummy_comp, num_oracle_offspring=5
    )

    # Select a mate for a given parent
    parent_idx = [0]  # Indices of the parents in the population
    # Mock the evaluation of each candidate
    selection.eval_mate = Mock(side_effect=[0.2, 0.8, 0.4])
    mate_idx = selection.select_mate(parent_idx, [3, 5, 7], problem, pop)

    # Check that the selected mate index is valid
    assert mate_idx == 1
    assert (
        selection.eval_mate.call_count == 3
    )  # Ensure eval_mate was called for each candidate
