from mating_kernel.pymoo.pref_model.oracle_pref_model import OraclePreferenceModel
from copy import deepcopy
from unittest.mock import Mock
from pymoo.problems import get_problem

from mating_kernel.problems.mo_tracking import make_tracked
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper
import pytest
import numpy as np
from mating_kernel.pymoo.mating.pref_tourn_sel import (
    PreferenceTournamentSelection,
)

from .. import dummy_comp, create_pop


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

    result = OraclePreferenceModel.count_survivors(
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

    mate_opt_idx = [2, 3, 4]  # Indices of the candidate mates in the population
    p_idx = [0]  # Indices of the parents in the population
    parents = [pop[i] for i in p_idx]  # Select parents based on indices
    mate_options = [pop[i] for i in mate_opt_idx]  # Select mate
    matings = [list(parents) + [mate] for mate in mate_options]
    # Create an instance of OraclePrefTournamentSelection
    selection = OraclePreferenceModel(func_comp=dummy_comp, num_oracle_simulations=5)

    # Simulate the matings with a fixed random seed
    survival_counts_1 = OraclePreferenceModel.simulate_matings(
        matings,
        problem,
        pop,
        selection.mating,
        selection.survival,
        random_state=np.random.default_rng(42),
    )

    # Simulate again with the same seed to check for reproducibility
    survival_counts_2 = OraclePreferenceModel.simulate_matings(
        matings,
        problem,
        pop,
        selection.mating,
        selection.survival,
        random_state=np.random.default_rng(42),
    )
    assert np.array_equal(survival_counts_1, survival_counts_2)


def test_simulate_matings():
    # Create a dummy problem and population
    problem = get_problem("zdt1")

    pop = create_pop(problem, size=10)
    mate_idx = [2, 3, 4]  # Indices of the candidate mates in the population
    p_idx = [0]  # Indices of the parents in the population
    parents = [pop[i] for i in p_idx]  # Select parents based on indices
    mate_options = [pop[i] for i in mate_idx]  # Select mate
    matings = [list(parents) + [mate] for mate in mate_options]
    # Create an instance of OraclePrefTournamentSelection
    selection = OraclePreferenceModel(func_comp=dummy_comp, num_oracle_simulations=5)

    # Simulate the matings
    survival_counts = OraclePreferenceModel.simulate_matings(
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


def test_count_survivors_per_mating_edge_cases():
    # Edge case: No offspring survived
    n_matings = 2
    n_offsprings = 2
    offspring = [object() for _ in range(4)]
    new_pop = []
    result = OraclePreferenceModel.count_survivors(
        offspring,
        new_pop,
        n_matings,
        n_offsprings,
    )
    assert result == [0, 0]

    # Edge case: All offspring survived
    new_pop = offspring.copy()
    result = OraclePreferenceModel.count_survivors(
        offspring,
        new_pop,
        n_matings,
        n_offsprings,
    )
    assert result == [2, 2]


def test_evaluate():
    num_sim = 5  # Number of simulations for mate selection
    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)
    pref_model = OraclePreferenceModel(
        func_comp=dummy_comp, num_oracle_simulations=num_sim
    )
    p_idx = [0]  # Indices of the parents in the population
    mate_idx = [2, 3, 4]  # Indices of the candidate mates in the population

    # Evaluate the mate options
    parents = [pop[i] for i in p_idx]  # Select parents based on indices
    mate_options = [pop[i] for i in mate_idx]  # Select mate
    scores = pref_model._evaluate(parents, mate_options, problem, pop)

    assert len(scores) == len(mate_options)
    assert all(isinstance(score, float) for score in scores)
    assert all(score >= pref_model.xl for score in scores)  # Check lower bound
    assert all(score <= pref_model.xu for score in scores)  # Check upper bound


@pytest.mark.parametrize("mock", [True, False])
def test_integration(mock):
    # Create a dummy problem and population
    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)
    num_sim = 5  # Number of simulations for mate selection

    pref_model = OraclePreferenceModel(
        func_comp=dummy_comp, num_oracle_simulations=num_sim
    )
    selection = PreferenceTournamentSelection(
        func_comp=dummy_comp,
        preference_model=pref_model,
        num_candidates=3,
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
        pref_model.simulate_matings = Mock(side_effect=side_effect)
    mate_idx = selection.select_mate(
        parent_idx, [3, 5, 7], problem, pop, random_state=None
    )

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
            pref_model.simulate_matings.call_count == num_sim
        )  # Ensure eval_mate was called for each candidate


def test_archive_disable():
    tracked_REProblem = make_tracked(REProblem)
    static_problem = tracked_REProblem(REProblemConfig(), inst=0)
    problem = PymooProblemWrapper(static_problem)
    pop = create_pop(problem, size=10)
    pop_cpy = deepcopy(pop)  # Make a copy of the original population
    num_sim = 5  # Number of simulations for mate selection
    len_records = len(problem.records)
    stats = problem.static_problem.get_archive_stats()
    pref_model = OraclePreferenceModel(
        func_comp=dummy_comp, num_oracle_simulations=num_sim
    )
    selection = PreferenceTournamentSelection(
        func_comp=dummy_comp,
        preference_model=pref_model,
        num_candidates=3,
    )
    selection._do(
        problem, pop, n_select=5, n_parents=2, random_state=np.random.default_rng(42)
    )
    # Check that the original population has not been modified
    assert len(pop) == 10
    for ind_original, ind_copied in zip(pop, pop_cpy):
        assert ind_original.F.tolist() == ind_copied.F.tolist()
        assert ind_original.data == ind_copied.data

    # Ensure that the archive and records have not been modified
    assert len(problem.records) == len_records
    assert problem.static_problem.get_archive_stats() == stats
    assert problem.record is True  # Ensure that the record attribute is still True
    assert (
        problem.static_problem._enabled is True
    )  # Ensure that the _enabled attribute is still True
