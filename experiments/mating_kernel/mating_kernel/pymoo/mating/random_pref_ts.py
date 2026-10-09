import numpy as np
from typing import Sequence

from mating_kernel.pymoo.mating.mating_pref_tournament_selection import (
    PreferredMatingTournamentSelection,
)
from pymoo.core.problem import Problem
from pymoo.core.population import Population


class RandomPrefTournamentSelection(PreferredMatingTournamentSelection):

    def select_mate(
        self,
        parents: Sequence[int],
        options: Sequence[int],
        problem: Problem,
        pop: Population,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ) -> int:
        # Select a random mate from the population
        if random_state is None:
            rng = np.random.default_rng()
            mate_idx = rng.integers(low=0, high=len(options))
        else:
            mate_idx = random_state.integers(low=0, high=len(options))
        return int(mate_idx)
