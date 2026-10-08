from pymoo.operators.selection.tournament import TournamentSelection
import numpy as np
from abc import ABC, abstractmethod


class PreferredMatingTournamentSelection(TournamentSelection, ABC):

    def _do(
        self,
        problem,
        pop,
        n_select: int,
        n_parents: int = 2,
        random_state: np.random.Generator | None = None,
        **kwargs
    ) -> np.ndarray:
        # Perform tournament selection to select individuals from the population
        selection = super()._do(
            problem, pop, n_select, n_parents, random_state=random_state, **kwargs
        )
        # Replace the last parent in each  mating with the preferred mate selected using the select_mate method
        for i in range(n_select):
            parents = selection[i].tolist()
            mate_idx = self.select_mate(
                parents, pop, problem, random_state=random_state, **kwargs
            )
            selection[i][-1] = mate_idx
        return selection

    @abstractmethod
    def select_mate(self, parents, pop, problem, random_state, **kwargs) -> int: ...
