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
import pytest
import numpy as np


@pytest.mark.parametrize("parents", [[0], [1, 2]])
def test_generate_matings(parents):
    # Create a dummy problem and population
    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)
    pop_cpy = deepcopy(pop)  # Make a copy of the original population

    mate_options = [2, 3, 4]  # Indices of the candidate mates in the population
    matings = OraclePrefTournamentSelection.generate_matings(parents, mate_options, pop)

    assert len(matings) == len(mate_options)
    assert all(
        pop != pop_cpy
    )  # Ensure the original population is not connected to these individuals anymore


@pytest.mark.parametrize(
    ("survival_flags", "expected"),
    [
        ([1, 0, 1, 1, 1, 0], [2, 1, 1]),
        ([0, 0, 0, 0, 0, 0], [0, 0, 0]),
        ([1, 1, 1, 1, 1, 1], [2, 2, 2]),
        ([1, 0, 0, 0, 1, 1], [1, 1, 1]),
    ],
)
def test_count_survivors_per_mating(survival_flags, expected):
    n_matings = 3
    n_offsprings = 2

    offspring = [object() for _ in range(6)]
    new_pop = [child for child, survived in zip(offspring, survival_flags) if survived]

    result = OraclePrefTournamentSelection.count_survivors(
        offspring,
        new_pop,
        n_matings,
        n_offsprings,
    )

    assert result == expected


def test_simulate_matings_seeding():
    # Create a dummy problem and population
    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)

    mate_options = [2, 3, 4]  # Indices of the candidate mates in the population
    parents = [0]  # Indices of the parents in the population
    matings = OraclePrefTournamentSelection.generate_matings(parents, mate_options, pop)

    # Create an instance of OraclePrefTournamentSelection
    selection = OraclePrefTournamentSelection(
        func_comp=dummy_comp, num_oracle_simulations=5
    )

    # Simulate the matings with a fixed random seed
    survival_counts_1 = OraclePrefTournamentSelection.simulate_matings(
        matings,
        problem,
        pop,
        selection.mating,
        selection.survival,
        random_state=np.random.default_rng(42),
    )

    # Simulate again with the same seed to check for reproducibility
    survival_counts_2 = OraclePrefTournamentSelection.simulate_matings(
        matings,
        problem,
        pop,
        selection.mating,
        selection.survival,
        random_state=np.random.default_rng(42),
    )

    assert survival_counts_1 == survival_counts_2


@pytest.mark.parametrize("archive_enabled", [True, False])
def test_simulate_matings(archive_enabled):
    # Create a dummy problem and population
    if not archive_enabled:
        problem = get_problem("zdt1")
    else:
        tracked_REProblem = make_tracked(REProblem)
        static_problem = tracked_REProblem(REProblemConfig(), inst=0)
        problem = PymooProblemWrapper(static_problem)
    pop = create_pop(problem, size=10)
    pop_cpy = deepcopy(pop)  # Make a copy of the original population
    if archive_enabled:
        stats = (
            problem.static_problem.get_archive_stats()
        )  # Ensure the archive is initialized
        assert stats["size"] > 0  # Previous evaluations
        len_records = len(problem.records)  # Store the original length of records

    mate_options = [2, 3, 4]  # Indices of the candidate mates in the population
    parents = [0]  # Indices of the parents in the population
    matings = OraclePrefTournamentSelection.generate_matings(parents, mate_options, pop)

    # Create an instance of OraclePrefTournamentSelection
    selection = OraclePrefTournamentSelection(
        func_comp=dummy_comp, num_oracle_simulations=5
    )

    # Simulate the matings
    survival_counts = OraclePrefTournamentSelection.simulate_matings(
        matings,
        problem,
        pop,
        selection.mating,
        selection.survival,
        random_state=np.random.default_rng(42),
    )

    assert len(survival_counts) == len(mate_options)
    assert all(
        count >= 0 for count in survival_counts
    )  # Ensure that survival counts are non-negative
    assert all(
        count <= selection.mating.crossover.n_offsprings for count in survival_counts
    )

    # Check that the original population has not been modified
    assert len(pop) == 10
    for ind_original, ind_copied in zip(pop, pop_cpy):
        assert ind_original.F.tolist() == ind_copied.F.tolist()
        assert ind_original.data == ind_copied.data

    if archive_enabled:
        # Ensure that the archive and records have not been modified
        assert len(problem.records) == len_records
        assert problem.static_problem.get_archive_stats() == stats


@pytest.mark.parametrize("mock", [True, False])
def test_select_mate(mock):
    # Create a dummy problem and population
    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)
    num_sim = 5  # Number of simulations for mate selection

    # Create an instance of OraclePrefTournamentSelection
    selection = OraclePrefTournamentSelection(
        func_comp=dummy_comp, num_oracle_simulations=num_sim
    )

    # Select a mate for a given parent
    parent_idx = [0]  # Indices of the parents in the population
    side_effect = [
        [1, 2, 0],
        [0, 1, 2],
        [2, 0, 1],
        [1, 1, 0],
        [0, 2, 1],
    ]

    if mock:
        # Mock the evaluation of each candidate
        selection.simulate_matings = Mock(side_effect=side_effect)
    mate_idx = selection.select_mate(parent_idx, [3, 5, 7], problem, pop)

    if not mock:
        # If not mocked, we can't predict the mate index, so we just check that it's valid
        assert 0 <= mate_idx < 3
    else:
        # Check that the selected mate index is valid
        assert mate_idx == 1
        side_effect_score_per_mate = [sum(x) for x in zip(*side_effect)]
        expected_selected_mate = side_effect_score_per_mate.index(
            max(side_effect_score_per_mate)
        )
        assert mate_idx == expected_selected_mate

        assert (
            selection.simulate_matings.call_count == num_sim
        )  # Ensure eval_mate was called for each candidate
