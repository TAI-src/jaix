import copy
from collections.abc import Callable, Sequence

import numpy as np
from pymoo.core.population import Population
from pymoo.core.problem import Problem
from pymoo.operators.selection.tournament import TournamentSelection

from mating_kernel.problems.mo_tracking import MOTrackingMixin
from mating_kernel.pymoo.pref_model.pref_model import PreferenceModel
from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper


class PreferenceTournamentSelection(TournamentSelection):

    def __init__(
        self,
        func_comp: Callable,
        preference_model: PreferenceModel,
        num_candidates: int = 1,
        candidate_pressure: int = 2,  # This is the default in pymoo
        **kwargs,
    ):
        super().__init__(func_comp=func_comp, **kwargs)
        self.num_candidates = num_candidates
        self.candidate_selection = TournamentSelection(
            func_comp=func_comp, pressure=candidate_pressure
        )
        self.preference_model = preference_model

    def _do(
        self,
        problem: Problem,
        pop: Population,
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

        # Make copies of pop and problem to prevent any modifications to the original objects
        prob_cpy = copy.deepcopy(problem)
        if isinstance(prob_cpy, PymooProblemWrapper):
            prob_cpy.record = False
            if isinstance(prob_cpy.static_problem, MOTrackingMixin):
                prob_cpy.static_problem.disable_adding()
        elif isinstance(prob_cpy, MOTrackingMixin):
            prob_cpy.disable_adding()
        pop_cpy = copy.deepcopy(pop)

        # Replace the last parent in each  mating with the preferred mate selected using the select_mate method
        for i in range(n_select):
            parents = selection[i].tolist()
            # Remove the last parent from the list of parents since that will be
            # replaced by the selected mate
            parents.pop(-1)
            options = preselect[i].tolist()
            mate_idx = self.select_mate(
                parents, options, prob_cpy, pop_cpy, random_state=random_state, **kwargs
            )
            selection[i][-1] = options[mate_idx]
        return selection

    def select_mate(
        self,
        parents: Sequence[int],
        options: Sequence[int],
        problem: Problem,
        pop: Population,
        random_state: np.random.Generator | None,
        **kwargs,
    ) -> int:
        """
        Returns the index of the selected mate from the options list based on the parents and problem context.
        This method should be implemented in subclasses to define the specific mate selection strategy.
        """
        scores = self.preference_model.evaluate(
            parents=[pop[i] for i in parents],
            mate_options=[pop[i] for i in options],
            problem=problem,
            pop=pop,
            random_state=random_state,
            **kwargs,
        )
        return int(np.argmax(scores))
