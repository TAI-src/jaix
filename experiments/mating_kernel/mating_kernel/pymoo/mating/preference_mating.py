import math
from copy import deepcopy

import numpy as np
from pymoo.core.crossover import Crossover
from pymoo.core.mating import Mating
from pymoo.core.mutation import Mutation
from pymoo.core.population import Population
from pymoo.core.problem import Problem
from pymoo.operators.selection.tournament import TournamentSelection

from mating_kernel.pymoo.mating.pref_tourn_sel import PreferenceTournamentSelection
from mating_kernel.pymoo.pref_model.pref_model import PreferenceModel


class PreferenceMating(Mating):
    def __init__(
        self,
        selection: TournamentSelection,
        crossover: Crossover,
        mutation: Mutation,
        pref_model: PreferenceModel,
        num_candidates: int = 1,
        candidate_pressure: int = 2,
        **kwargs,
    ):
        assert isinstance(
            selection, TournamentSelection
        ), "selection must be an instance of TournamentSelection"
        super().__init__(selection, crossover, mutation, **kwargs)
        self.pref_model = pref_model
        self.selection = PreferenceTournamentSelection(
            func_comp=selection.func_comp,
            preference_model=pref_model,
            num_candidates=num_candidates,
            candidate_pressure=candidate_pressure,
        )

    def _set_pref(
        self,
        parents: Population,
        off: Population,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ):
        # Create a new population of parents with their preference vectors as their decision variables for the preference model

        pref_parents = deepcopy(parents)
        for mating in pref_parents:
            for p in mating:
                p.X = p.pref

        # The dummy problem is just used for information about the preference model, such as the number of preference variables and their bounds.
        # It is not used for any optimization or evaluation.
        off_pref = self.crossover(
            self.pref_model.dummy_problem,
            pref_parents,
            random_state=random_state,
            **kwargs,
        )
        off_pref = self.mutation(
            self.pref_model.dummy_problem, off_pref, random_state=random_state, **kwargs
        )

        # Now assign the preferences of the offspring created through crossover and mutation to the offspring created through crossover and mutation
        for i in range(len(off)):
            off[i].pref = off_pref[i].X

    @staticmethod
    def init_pref(
        pref_model: PreferenceModel,
        pop: Population,
        random_state: np.random.Generator | None = None,
    ):
        """
        Init the preference vector if not already exists"""
        if not pref_model.pref_inheritance:
            return  # If preference inheritance is disabled, do not initialize preferences
        for ind in pop:
            pref_model.init_individual(ind, random_state=random_state)

    def create_offspring(
        self,
        problem: Problem,
        parents: Population,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ) -> Population:
        off = self.crossover(problem, parents, random_state=random_state, **kwargs)
        # do the mutation on the offsprings created through crossover
        off = self.mutation(problem, off, random_state=random_state, **kwargs)
        return off

    def _do(
        self,
        problem: Problem,
        pop: Population,
        n_offsprings: int,
        parents: np.ndarray | None = None,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ):
        # how many parents need to be select for the mating - depending on number of offsprings remaining
        n_matings = math.ceil(n_offsprings / self.crossover.n_offsprings)
        PreferenceMating.init_pref(self.pref_model, pop, random_state=random_state)

        # if the parents for the mating are not provided directly - usually selection will be used
        if parents is None:
            # select the parents for the mating - just an index array
            parent_pop: Population = self.selection(
                problem,
                pop,
                n_matings,
                n_parents=self.crossover.n_parents,
                random_state=random_state,
                **kwargs,
            )
        # Apply crossover and mutation to the selected parents to create the offspring
        # This is like the regular mating, operating on X
        off = self.create_offspring(
            problem, parent_pop, random_state=random_state, **kwargs
        )

        if self.pref_model.pref_inheritance:
            # If preference inheritance is enabled, we will perform crossover and mutation on the preference vectors of the parents to create the preference vectors for the offspring
            self._set_pref(parent_pop, off, random_state=random_state, **kwargs)
        return off
