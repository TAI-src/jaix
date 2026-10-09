from pymoo.core.mating import Mating
from pymoo.operators.selection.tournament import TournamentSelection
import math
from copy import deepcopy

from mating_kernel.pymoo.pref_model.pref_model import PreferenceModel
from mating_kernel.pymoo.mating.pref_tourn_sel import PreferenceTournamentSelection


class PreferenceMating(Mating):
    def __init__(
        self,
        selection,
        crossover,
        mutation,
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

    def _do_pref(self, parents, off, random_state=None, **kwargs):
        # TODO: This is just a stub for now until other preference models are implemented.
        # So this is never called
        pref_parents = []
        for mating in parents:
            mates = [deepcopy(parent) for parent in mating]
            for p in mates:
                p.X = p.data["preferences"]
            pref_parents.append(mates)

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
            off[i].data["preferences"] = off_pref[i].X

    def _do(
        self, problem, pop, n_offsprings, parents=None, random_state=None, **kwargs
    ):
        # how many parents need to be select for the mating - depending on number of offsprings remaining
        n_matings = math.ceil(n_offsprings / self.crossover.n_offsprings)

        # if the parents for the mating are not provided directly - usually selection will be used
        if parents is None:
            # select the parents for the mating - just an index array
            parents = self.selection(
                problem,
                pop,
                n_matings,
                n_parents=self.crossover.n_parents,
                random_state=random_state,
                **kwargs,
            )
        # Apply crossover and mutation to the selected parents to create the offspring
        # This is like the regular mating, operating on X
        off = self.crossover(problem, parents, random_state=random_state, **kwargs)
        # do the mutation on the offsprings created through crossover
        off = self.mutation(problem, off, random_state=random_state, **kwargs)

        return off
