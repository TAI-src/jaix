from abc import ABC, abstractmethod

import numpy as np
from pymoo.operators.selection.tournament import TournamentSelection


class PreferredMatingTournamentSelection(TournamentSelection, ABC):

    def __init__(
        self,
        func_comp,
        num_candidates: int = 1,
        candidate_pressure: int = 2,  # This is the default in pymoo
        **kwargs,
    ):
        super().__init__(func_comp=func_comp, **kwargs)
        self.num_candidates = num_candidates
        self.candidate_selection = TournamentSelection(
            func_comp=func_comp, pressure=candidate_pressure
        )

    def _do(
        self,
        problem,
        pop,
        n_select: int,
        n_parents: int = 2,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ) -> np.ndarray:
        # Perform tournament selection to select individuals from the population
        selection = super()._do(
            problem, pop, n_select, n_parents, random_state=random_state, **kwargs
        )
        preselect = self.candidate_selection._do(
            problem,
            pop,
            n_select,
            self.num_candidates,
            random_state=random_state,
            **kwargs,
        )

        # Replace the last parent in each  mating with the preferred mate selected using the select_mate method
        for i in range(n_select):
            parents = selection[i].tolist()
            options = preselect[i].tolist()
            mate_idx = self.select_mate(
                parents, options, problem, random_state=random_state, **kwargs
            )
            selection[i][-1] = options[mate_idx]
        return selection

    @abstractmethod
    def select_mate(self, parents, options, problem, random_state, **kwargs) -> int: ...
