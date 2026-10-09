from mating_kernel.pymoo.mating.pref_tourn_sel import PreferenceTournamentSelection
import numpy as np
from pymoo.operators.selection.tournament import TournamentSelection
from mating_kernel.pymoo.pref_model.random_pref_model import RandomPreferenceModel

from .. import create_pop, dummy_comp


def test_tournament_selection_with_random_pref():
    # Create a RandomPrefTournamentSelection instance
    original_selection = TournamentSelection(func_comp=dummy_comp)
    selection = PreferenceTournamentSelection(
        func_comp=dummy_comp,
        preference_model=RandomPreferenceModel(),
        num_candidates=3,
        candidate_pressure=5,
    )

    # Create a population of individuals
    pop = create_pop(None, size=10)
    off2 = original_selection.do(
        None, pop, n_select=3, n_parents=2, random_state=np.random.default_rng(42)
    )
    off = selection.do(
        None, pop, n_select=3, n_parents=2, random_state=np.random.default_rng(42)
    )

    assert off.shape == (3, 2)  # Check that the output shape is correct
    # Check that the first individuals in each mating are the same as in the original selection
    # and that the last one is not necessarily
    different = 0
    for orig_mates, new_mates in zip(off2, off):
        assert orig_mates[0] == new_mates[0]
        if orig_mates[1] != new_mates[1]:
            different += 1
    assert different > 0  # Ensure that at least one mate is different
