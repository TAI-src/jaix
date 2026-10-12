from mating_kernel.pymoo.pref_model.pref_model import PreferenceModel
import numpy as np
from typing import Sequence
import pytest
from pymoo.core.individual import Individual


class DummyPreferenceModel(PreferenceModel):

    def __init__(self, pref_inheritance: bool = True, n_preferences: int = 1) -> None:
        self.pref_inheritance = pref_inheritance
        self._n_preferences = n_preferences
        self.calls = []  # type: ignore

    @property
    def n_preferences(self) -> int:
        return self._n_preferences

    @property
    def xl(self) -> np.ndarray:
        return np.array([0.0] * self._n_preferences)

    @property
    def xu(self) -> np.ndarray:
        return np.array([1.0] * self._n_preferences)

    def _evaluate(
        self,
        parents,
        mate_options,
        problem,
        pop,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ) -> Sequence[float]:
        self.calls.append((parents, mate_options))
        vals = np.zeros(len(mate_options))
        vals[1] = 1.0  # Always prefer the second option
        return list(vals)


def test_dummy_preference_model():
    model = DummyPreferenceModel()
    parents = [0, 1]
    mate_options = [2, 3, 4]
    problem = None
    pop = None
    random_state = np.random.default_rng(42)

    # Call the _evaluate method
    preferences = model._evaluate(parents, mate_options, problem, pop, random_state)

    # Check that the preferences are as expected
    assert len(preferences) == len(mate_options)
    assert preferences[0] == 0.0
    assert preferences[1] == 1.0
    assert preferences[2] == 0.0

    # Check that the calls list has been updated correctly
    assert len(model.calls) == 1
    assert model.calls[0][0] == parents
    assert model.calls[0][1] == mate_options


@pytest.mark.parametrize("n_preferences", [1, 2, 5])
def test_initialize(n_preferences):
    pref_model = DummyPreferenceModel(n_preferences=n_preferences)
    init = pref_model.initialize(np.random.default_rng(42))
    assert len(init) == pref_model.n_preferences


@pytest.mark.parametrize("pref_inheritance", [True, False])
@pytest.mark.parametrize("n_preferences", [1, 2, 5])
def test_init_individual(pref_inheritance, n_preferences):
    pref_model = DummyPreferenceModel(
        pref_inheritance=pref_inheritance, n_preferences=n_preferences
    )
    ind = Individual()
    pref_model.init_individual(ind, random_state=np.random.default_rng(42))

    if pref_inheritance:
        assert hasattr(ind, "pref")
        assert len(ind.pref) == pref_model.n_preferences
    else:
        assert not hasattr(ind, "pref") or ind.pref is None


def test_init_individual_overwrite():
    pref_model = DummyPreferenceModel(pref_inheritance=True, n_preferences=3)
    ind = Individual()
    ind.pref = np.array([0.5, 0.5, 0.5])  # Pre-existing preference vector
    pref_model.init_individual(ind, random_state=np.random.default_rng(42))

    # The preference vector should not be overwritten
    assert np.array_equal(ind.pref, np.array([0.5, 0.5, 0.5]))


@pytest.mark.parametrize("n_preferences", [1, 2, 5])
def test_problem(n_preferences):
    pref_model = DummyPreferenceModel(n_preferences=n_preferences)
    problem = pref_model.dummy_problem

    assert problem.n_var == pref_model.n_preferences
    assert problem.n_obj == 1
    assert problem.n_constr == 0
    assert len(problem.xl) == pref_model.n_preferences
    assert len(problem.xu) == pref_model.n_preferences
