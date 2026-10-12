from mating_kernel.pymoo.mating.preference_mating import PreferenceMating
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.selection.tournament import TournamentSelection
from pymoo.problems import get_problem
from pymoo.algorithms.moo.nsga2 import NSGA2, binary_tournament
from pymoo.optimize import minimize
from ..pref_model.test_pref_model import DummyPreferenceModel
import pytest

from .. import create_pop, dummy_comp
from pymoo.core.individual import Individual
import numpy as np
import copy


def test_do_without_prefs():
    selection = TournamentSelection(func_comp=dummy_comp)
    crossover = SBX(prob=0.9, eta=15)
    mutation = PM(eta=20)
    mating = PreferenceMating(
        selection=selection,
        crossover=crossover,
        mutation=mutation,
        pref_model=DummyPreferenceModel(pref_inheritance=False),
        num_candidates=5,
    )

    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)
    off = mating.do(problem, pop, n_offsprings=5, random_state=None)

    assert len(off) == 5
    assert all(isinstance(ind, Individual) for ind in off)


@pytest.mark.parametrize("pref_inheritance", [True, False])
def test_init_pref(pref_inheritance):
    pref_model = DummyPreferenceModel(
        pref_inheritance=pref_inheritance, n_preferences=3
    )
    pop = create_pop(None, size=10)
    PreferenceMating.init_pref(pref_model, pop, random_state=np.random.default_rng(42))
    for ind in pop:
        if pref_inheritance:
            assert hasattr(ind, "pref")
            assert len(ind.pref) == pref_model.n_preferences
        else:
            assert not hasattr(ind, "pref") or ind.pref is None


@pytest.mark.parametrize("n_preferences", [0, 1, 5])
def test_create_offspring(n_preferences):
    selection = TournamentSelection(func_comp=dummy_comp)
    crossover = SBX(prob=0.9, eta=15)
    mutation = PM(eta=20)
    pref_model = DummyPreferenceModel(
        pref_inheritance=True, n_preferences=n_preferences
    )
    mating = PreferenceMating(
        selection=selection,
        crossover=crossover,
        mutation=mutation,
        pref_model=pref_model,
        num_candidates=5,
    )

    if n_preferences > 0:
        problem = pref_model.dummy_problem
    else:
        problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)
    # Create 3 matings of 2 parents each
    parent_idx = np.random.default_rng(42).choice(len(pop), size=(3, 2), replace=False)
    parent_pop = [[pop[i] for i in idx_pair] for idx_pair in parent_idx]

    off = mating.create_offspring(
        problem, parent_pop, random_state=np.random.default_rng(42)
    )

    assert len(off) == len(parent_pop) * crossover.n_offsprings
    assert all(isinstance(ind, Individual) for ind in off)


@pytest.mark.parametrize("n_preferences", [1, 2, 5])
def test_do_with_prefs(n_preferences):
    selection = TournamentSelection(func_comp=dummy_comp)
    crossover = SBX(prob=0.9, eta=15)
    mutation = PM(eta=20)
    pref_model = DummyPreferenceModel(
        pref_inheritance=True, n_preferences=n_preferences
    )
    mating = PreferenceMating(
        selection=selection,
        crossover=crossover,
        mutation=mutation,
        pref_model=pref_model,
        num_candidates=5,
    )

    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)  # This represent the original population
    PreferenceMating.init_pref(pref_model, pop, random_state=np.random.default_rng(42))
    pop_copy = copy.deepcopy(pop)  # Store a copy of the original population
    off = create_pop(
        problem, size=6
    )  # These simulate the offspring of the regular variation
    x_copy_off = [ind.X.copy() for ind in off]  # Store original decision variables
    off_copy = copy.deepcopy(off)  # Store a copy of the original offspring
    # Create 3 matings of 2 parents each
    parent_idx = np.random.default_rng(42).choice(len(pop), size=(3, 2), replace=False)
    parent_pop = [[pop[i] for i in idx_pair] for idx_pair in parent_idx]

    mating._set_pref(parent_pop, off, random_state=np.random.default_rng(42))
    for x_copy, orig_off, ind in zip(x_copy_off, off_copy, off):
        assert np.array_equal(x_copy, ind.X)  # Ensure decision variables are unchanged
        assert np.array_equal(
            orig_off.X, ind.X
        )  # Ensure decision variables are unchanged
        assert np.array_equal(
            orig_off.F, ind.F
        )  # Ensure objective values are unchanged
        assert orig_off.data == ind.data  # Ensure data is unchanged
    # Ensure that the preferences have been set correctly
    for ind in off:
        assert hasattr(ind, "pref")
        assert len(ind.pref) == n_preferences
    # Ensure that the preferences of the parents have not been modified
    for orig_ind, now_ind in zip(pop, pop_copy):
        assert np.array_equal(orig_ind.X, now_ind.X)
        assert np.array_equal(orig_ind.F, now_ind.F)
        assert orig_ind.data == now_ind.data
        assert np.array_equal(orig_ind.pref, now_ind.pref)


@pytest.mark.parametrize("pref_inheritance", [True, False])
@pytest.mark.parametrize("n_preferences", [1, 2, 5])
def test_integration_with_nsga2(pref_inheritance, n_preferences):
    selection = TournamentSelection(func_comp=binary_tournament)
    crossover = SBX(prob=0.9, eta=15)
    mutation = PM(eta=20)
    mating = PreferenceMating(
        selection=selection,
        crossover=crossover,
        mutation=mutation,
        pref_model=DummyPreferenceModel(
            pref_inheritance=pref_inheritance, n_preferences=n_preferences
        ),
        num_candidates=5,
    )

    problem = get_problem("zdt1")
    algorithm = NSGA2(
        mating=mating,
    )
    res = minimize(
        problem,
        algorithm,
        termination=("n_gen", 5),
        seed=42,
        save_history=True,
    )

    assert res.history is not None
    assert len(res.history) == 5  # 5 generations
