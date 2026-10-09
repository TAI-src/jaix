from mating_kernel.pymoo.pref_model.pref_model import PreferenceModel
from pymoo.core.individual import Individual
from typing import Sequence
import numpy as np
from pymoo.core.problem import Problem
from pymoo.core.population import Population


class RandomPreferenceModel(PreferenceModel):
    pref_inheritance = False  # Preferences are not inherited from parents

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
        parents: Sequence[Individual],
        mate_options: Sequence[Individual],
        problem: Problem,
        pop: Population,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ) -> Sequence[float]:
        # Return random preference scores for each mate option
        scores = (
            random_state.uniform(0, 1, size=len(mate_options))
            if random_state
            else np.random.uniform(0, 1, size=len(mate_options))
        )
        return list(scores)
