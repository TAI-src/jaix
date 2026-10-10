import copy
from collections.abc import Callable, Sequence

import numpy as np
from mating_kernel.problems.mo_tracking import MOTrackingMixin
from mating_kernel.pymoo.pref_model.pref_model import PreferenceModel
from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper
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


class OraclePreferenceModel(PreferenceModel):
    pref_inheritance = False  # Preferences are not inherited from parents

    def __init__(
        self,
        func_comp: Callable,
        num_oracle_simulations: int = 30,
        crossover=None,
        mutation=None,
        survival=None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        selection = TournamentSelection(func_comp=func_comp)

        crossover = crossover or SBX(prob=0.9, eta=15)  # Default from pymoo nsga2
        mutation = mutation or PM(eta=20)  # Default from pymoo nsga2
        survival = survival or RankAndCrowdingSurvival()  # Default from pymoo nsga2
        self.mating = Mating(
            selection=selection, crossover=crossover, mutation=mutation
        )
        self.survival = survival
        assert (
            num_oracle_simulations > 0
        ), "num_oracle_simulations must be greater than 0"
        self.num_oracle_simulations = num_oracle_simulations

    @property
    def n_preferences(self) -> int:
        return 1

    @property
    def xl(self) -> np.ndarray:
        return np.array([0.0])

    @property
    def xu(self) -> np.ndarray:
        return np.array([1.0])

    @staticmethod
    def generate_matings(
        parents: Sequence[Individual],
        mate_options: Sequence[Individual],
        pop: Population,
    ) -> tuple[Sequence[Sequence[Individual]], Population]:
        """
        Creates a list of matings for all candidate mates.
        Individuals belong to a copied population to avoid modifying the original.
        """
        pop_cpy = copy.deepcopy(pop)

        # Map original individual identities to their population indices.
        index_by_id = {id(ind): i for i, ind in enumerate(pop)}

        parent_idx = [index_by_id[id(ind)] for ind in parents]

        matings = []
        for mate in mate_options:
            p = [pop_cpy[i] for i in parent_idx]
            mate_idx = index_by_id[id(mate)]
            p.append(pop_cpy[mate_idx])
            matings.append(p)

        return matings, pop_cpy

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
    ) -> Sequence[float]:

        off = mating._do(
            problem,
            pop,
            n_offsprings=-1,  # This is ignored
            parents=np.array(matings),
            random_state=random_state,
        )
        assert (
            len(off) == len(matings) * mating.crossover.n_offsprings
        ), "Unexpected number of offspring generated."
        Evaluator().eval(problem, off)

        new_pop = survival.do(
            problem,
            Population.merge(pop, off),
            n_survive=len(pop),
            random_state=random_state,
        )
        # Determine which offspring survived for each candidate mate
        survived_per_mate = OraclePreferenceModel.count_survivors(
            off,
            new_pop,
            n_matings=len(matings),
            n_offsprings=mating.crossover.n_offsprings,
        )
        survival_rate = np.array(survived_per_mate) / mating.crossover.n_offsprings

        return survival_rate

    def _evaluate(
        self,
        parents: Sequence[Individual],
        mate_options: Sequence[Individual],
        problem: Problem,
        pop: Population,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ) -> Sequence[float]:
        # Generate matings for all candidate mates
        matings, pop_cpy = self.generate_matings(parents, mate_options, pop)
        prob_cpy = copy.deepcopy(problem)
        if isinstance(prob_cpy, PymooProblemWrapper):
            prob_cpy.record = False
            if isinstance(prob_cpy.static_problem, MOTrackingMixin):
                prob_cpy.static_problem.disable_adding()
        elif isinstance(prob_cpy, MOTrackingMixin):
            prob_cpy.disable_adding()
        mate_scores = np.zeros(len(mate_options))
        for _ in range(self.num_oracle_simulations):
            survival_rates = self.simulate_matings(
                matings,
                prob_cpy,
                pop_cpy,
                self.mating,
                self.survival,
                random_state=random_state,
            )
            mate_scores += survival_rates
        # Average the survival rates over the number of simulations
        survival_rates = list(mate_scores / self.num_oracle_simulations)

        return survival_rates
