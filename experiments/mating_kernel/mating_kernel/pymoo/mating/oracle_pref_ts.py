import copy
from collections.abc import Callable, Sequence

import numpy as np
from pymoo.algorithms.moo.nsga2 import RankAndCrowdingSurvival
from pymoo.core.evaluator import Evaluator
from pymoo.core.individual import Individual
from pymoo.core.mating import Mating
from pymoo.core.population import Population
from pymoo.core.problem import Problem
from pymoo.core.survival import Survival
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.selection.tournament import TournamentSelection

from mating_kernel.pymoo.mating.mating_pref_tournament_selection import (
    PreferredMatingTournamentSelection,
)


class OraclePrefTournamentSelection(PreferredMatingTournamentSelection):
    def __init__(
        self,
        func_comp: Callable,
        num_oracle_simulations: int = 30,
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
        assert (
            num_oracle_simulations > 0
        ), "num_oracle_simulations must be greater than 0"
        self.num_oracle_simulations = num_oracle_simulations

    @staticmethod
    def generate_matings(
        parents: Sequence[int],
        mate_options: Sequence[int],
        pop: Population,
    ) -> Sequence[Sequence[Individual]]:
        # Generate matings for all candidate mates
        matings = []
        for mate in mate_options:
            p = [pop[i] for i in parents]
            p.append(pop[mate])
            matings.append(p)
        return matings

    @staticmethod
    def count_survivors(
        offspring: Sequence[Individual],
        new_pop: Population,
        n_matings: int,
        n_offsprings: int,
    ) -> Sequence[int]:
        # Determine which offspring survived for each candidate mate
        all_survived = np.array([child in new_pop for child in offspring], dtype=int)
        survived_per_mate = all_survived.reshape(n_offsprings, n_matings).sum(axis=0)
        return survived_per_mate.tolist()

    @staticmethod
    def simulate_matings(
        matings: Sequence[Sequence[Individual]],
        problem: Problem,
        pop: Population,
        mating: Mating,
        survival: Survival,
        random_state=None,
    ) -> Sequence[int]:

        pop_cpy = copy.deepcopy(pop)
        prob_cpy = copy.deepcopy(problem)
        off = mating._do(
            prob_cpy,
            pop_cpy,
            n_offsprings=-1,  # This is ignored
            parents=np.array(matings),
            random_state=random_state,
        )
        assert (
            len(off) == len(matings) * mating.crossover.n_offsprings
        ), "Unexpected number of offspring generated."
        Evaluator().eval(prob_cpy, off)

        new_pop = survival.do(
            prob_cpy, Population.merge(pop_cpy, off), n_survive=len(pop_cpy)
        )
        # Determine which offspring survived for each candidate mate
        survived_per_mate = OraclePrefTournamentSelection.count_survivors(
            off,
            new_pop,
            n_matings=len(matings),
            n_offsprings=mating.crossover.n_offsprings,
        )

        return survived_per_mate

    def select_mate(
        self,
        parents: Sequence[int],
        options: Sequence[int],
        problem: Problem,
        pop: Population,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ) -> int:

        # Generate offspring for all matings and evaluate them
        pop_cpy = copy.deepcopy(pop)
        matings = self.generate_matings(parents, options, pop_cpy)

        mate_scores = [0] * len(options)
        for _ in range(self.num_oracle_simulations):
            survived = self.simulate_matings(
                matings, problem, pop_cpy, self.mating, self.survival, random_state
            )
            mate_scores = [score + s for score, s in zip(mate_scores, survived)]

        # Select the candidate with the highest score
        best_candidate_idx = int(np.argmax(mate_scores))
        return best_candidate_idx
