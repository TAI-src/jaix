from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.selection.tournament import TournamentSelection

from mating_kernel.pymoo.do_recorder import DoRecorderMixin, RecordingConfig
from mating_kernel.pymoo.recordable_object import make_recordable, make_recorded

rec_config = RecordingConfig(save_output=True, meta_fields=["pressure"])


def test_make_recorded():
    # Create a recorded version of NSGA2
    selection = TournamentSelection(func_comp=lambda a, b: a < b)
    recorded_selection = make_recorded(selection, recording_config=rec_config)
    assert isinstance(recorded_selection, DoRecorderMixin)
    assert issubclass(type(recorded_selection), DoRecorderMixin)


def test_make_recordable():
    # Create a recordable version of NSGA2
    RecordableNSGA2 = make_recordable(NSGA2)

    # Instantiate the recordable NSGA2 with 'selection' to be recorded
    algorithm = RecordableNSGA2(pop_size=92, record_args={"selection": rec_config})

    # Check that the selection operator is now a subclass of DoRecorderMixin
    assert isinstance(algorithm.mating.selection, DoRecorderMixin)
    assert issubclass(type(algorithm.mating.selection), DoRecorderMixin)
