from mating_kernel.pymoo.mating.pref_tourn_sel import PreferenceTournamentSelection
import numpy as np

from ..pref_model.test_pref_model import DummyPreferenceModel
from .. import dummy_comp, create_pop


def test_preferred_mate_selection():
    pop = create_pop(size=10)
    pref_model = DummyPreferenceModel()
    selector = PreferenceTournamentSelection(
        func_comp=dummy_comp,
        preference_model=pref_model,
        num_candidates=3,
        candidate_pressure=5,
    )

    result = selector._do(
        None,
        pop,
        n_select=10,
        n_parents=2,
        random_state=np.random.default_rng(42),
    )

    assert len(pref_model.calls) == 10

    for i, (parents, options) in enumerate(pref_model.calls):
        assert len(options) == 3
        expected_option = options[
            1
        ]  # Since DummyPreferenceModel always returns index 1
        chosen_index = result[i, -1]
        assert (
            pop[chosen_index] == expected_option
        )  # Check that the chosen mate is the expected one


def test_candidate_selection_configuration():
    pref_model = DummyPreferenceModel()
    selector = PreferenceTournamentSelection(
        func_comp=dummy_comp,
        preference_model=pref_model,
        num_candidates=3,
        candidate_pressure=5,
    )

    assert selector.num_candidates == 3
    assert selector.candidate_selection.pressure == 5
