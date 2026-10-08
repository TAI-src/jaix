from mating_kernel.pymoo.mating.mating_pref_tournament_selection import (
    PreferredMatingTournamentSelection,
)
import numpy as np


class RandomPrefTournamentSelection(PreferredMatingTournamentSelection):
    def select_mate(
        self,
        parents,
        pop,
        problem,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ) -> int:
        # Select a random mate from the population
        if random_state is None:
            rng = np.random.default_rng()
            mate_idx = rng.integers(low=0, high=len(pop))
        else:
            mate_idx = random_state.integers(low=0, high=len(pop))
        return int(mate_idx)
