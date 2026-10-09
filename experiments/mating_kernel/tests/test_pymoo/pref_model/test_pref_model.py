from mating_kernel.pymoo.pref_model.pref_model import PreferenceModel
import numpy as np
from typing import Sequence

from .. import create_pop


class DummyPreferenceModel(PreferenceModel):
    pref_inheritance = True

    def __init__(self) -> None:
        self.calls = []

    @property
    def n_preferences(self) -> int:
        return 1

    @property
    def xl(self) -> np.ndarray:
        return np.array([0.0])

    @property
    def xu(self) -> np.ndarray:
        return np.array([1.0])

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
        return vals


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


def test_preference_model_initialization():
    model = DummyPreferenceModel()
    pop = create_pop(size=5)
    model.evaluate(pop, pop, None, pop)  # Call evaluate to populate calls
    assert pop[0].pref is not None
    assert pop[0].pref.shape == (model.n_preferences,)
    assert len(model.calls) == 1
