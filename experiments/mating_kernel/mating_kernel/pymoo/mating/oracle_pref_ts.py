from typing import Sequence, Callable
import numpy as np
import copy

from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.core.mating import Mating
from pymoo.algorithms.moo.nsga2 import RankAndCrowdingSurvival
from pymoo.operators.selection.tournament import TournamentSelection
from pymoo.core.evaluator import Evaluator
from pymoo.core.population import Population
from pymoo.core.problem import Problem

from mating_kernel.pymoo.mating.mating_pref_tournament_selection import (
    PreferredMatingTournamentSelection,
)


class OraclePrefTournamentSelection(PreferredMatingTournamentSelection):
    def __init__(
        self,
        func_comp: Callable,
        num_oracle_offspring: int,
        num_candidates: int = 1,
        candidate_pressure: int = 2,  # This is the default in pymoo
        **kwargs,
    ):
        super().__init__(
            func_comp=func_comp,
            num_candidates=num_candidates,
            candidate_pressure=candidate_pressure,
            **kwargs,
        )
        # FIXME: Should adapt to the algorithms's operators instead of hardcoding them here
        crossover = SBX(eta=15, prob=0.9)  # Default from pymoo nsga2
        mutation = PM(eta=20)  # Default from pymoo nsga2
        selection = TournamentSelection(func_comp=func_comp)

        self.mating = Mating(
            selection=selection, crossover=crossover, mutation=mutation
        )
        self.survival = RankAndCrowdingSurvival()  # Default from pymoo nsga2
        self.num_offspring = num_oracle_offspring

    def eval_mate(
        self,
        parents: Sequence[int],
        mate: int,
        problem: Problem,
        pop: Population,
        random_state=None,
    ) -> float:
        # Generate offspring with the given mate and evaluate them
        pop_cpy = copy.deepcopy(pop)
        assert len(parents) == 1, "This method currently supports only one parent."
        p = [pop_cpy[i] for i in parents]
        p.append(pop_cpy[mate])
        off = self.mating.do(
            problem,
            pop_cpy,
            n_offsprings=self.num_offspring,
            parents=np.array([p]),
            random_state=random_state,
        )
        Evaluator().eval(problem, off)
        new_pop = self.survival.do(
            problem, Population.merge(pop_cpy, off), n_survive=len(pop_cpy)
        )
        # Count the surviving offspring for the given mate
        survivors = [o for o in off if o in new_pop]

        return (
            len(survivors) / self.num_offspring
        )  # Return the proportion of surviving offspring

    def select_mate(
        self,
        parents: Sequence[int],
        options: Sequence[int],
        problem: Problem,
        pop: Population,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ) -> int:
        mate_scores = []
        # Generate offspring with each candidate and evaluate them
        for mate in options:
            score = self.eval_mate(
                parents,
                mate,
                problem,
                pop,
                random_state=random_state,
            )
            mate_scores.append(score)
        # Select the candidate with the highest score
        best_candidate_idx = int(np.argmax(mate_scores))
        return best_candidate_idx
