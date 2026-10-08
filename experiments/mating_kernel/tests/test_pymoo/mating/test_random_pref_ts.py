import numpy as np
from pymoo.operators.selection.tournament import TournamentSelection
from pymoo.problems import get_problem
from mating_kernel.pymoo.mating.random_pref_ts import RandomPrefTournamentSelection
from pymoo.core.population import Population


def dummy_comp(pop, P, **kwargs):
    return P[:, 0]


def create_pop(problem, size=10):
    x = np.random.uniform(low=problem.xl, high=problem.xu, size=(size, len(problem.xl)))
    pop = Population.new("X", x)
    pop.set("F", problem.evaluate(pop.get("X")))
    return pop


def test_random_pref_ts():
    # Create a RandomPrefTournamentSelection instance
    selection = RandomPrefTournamentSelection(func_comp=dummy_comp)
    original_selection = TournamentSelection(func_comp=dummy_comp)

    # Create a population of individuals
    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)
    off2 = original_selection.do(
        None, pop, n_select=3, n_parents=2, random_state=np.random.default_rng(42)
    )
    off = selection.do(
        None, pop, n_select=3, n_parents=2, random_state=np.random.default_rng(42)
    )

    assert off.shape == (3, 2)  # Check that the output shape is correct
    # Check that the first individuals in each mating are the same as in the original selection
    # and that the last one is not necessarily
    different = 0
    for orig_mates, new_mates in zip(off2, off):
        assert orig_mates[0] == new_mates[0]
        if orig_mates[1] != new_mates[1]:
            different += 1
    assert different > 0  # Ensure that at least one mate is different


def test_random_pref_ts_with_preselect():
    # Test with more preselect options
    selection = RandomPrefTournamentSelection(func_comp=dummy_comp, num_candidates=3)
    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)
    off = selection.do(
        None, pop, n_select=3, n_parents=2, random_state=np.random.default_rng(42)
    )
    assert off.shape == (3, 2)  # Check that the output shape is correct
