from mating_kernel.pymoo.mating.mating_pref_tournament_selection import (
    PreferredMatingTournamentSelection,
)
import numpy as np
from .test_random_pref_ts import dummy_comp, create_pop
from pymoo.problems import get_problem


def test_preferred_mate_selection():
    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)
    calls = []

    class TestSelection(PreferredMatingTournamentSelection):
        def select_mate(self, parents, options, problem, pop, **kwargs):
            calls.append((parents, options))
            return 1

    selector = TestSelection(
        num_candidates=3, candidate_pressure=5, func_comp=dummy_comp
    )

    result = selector._do(
        problem,
        pop,
        n_select=10,
        n_parents=2,
        random_state=np.random.default_rng(42),
    )

    assert len(calls) == 10

    for i, (parents, options) in enumerate(calls):
        assert len(options) == 3
        assert result[i, -1] == options[1]


def test_candidate_selection_configuration():
    class TestSelection(PreferredMatingTournamentSelection):
        def select_mate(self, parents, options, problem, pop, **kwargs):
            return 1

    selector = TestSelection(
        num_candidates=3,
        candidate_pressure=5,
        func_comp=dummy_comp,
    )

    assert selector.num_candidates == 3
    assert selector.candidate_selection.pressure == 5
