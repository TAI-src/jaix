from mating_kernel.pymoo.recordable_object import make_recordable, make_recorded
from pymoo.algorithms.moo.nsga2 import NSGA2
from mating_kernel.pymoo.do_recorder import DoRecorderMixin
from pymoo.operators.selection.tournament import TournamentSelection


def test_make_recorded():
    # Create a recorded version of NSGA2
    selection = TournamentSelection(func_comp=lambda a, b: a < b)
    recorded_selection = make_recorded(selection)
    assert isinstance(recorded_selection, DoRecorderMixin)
    assert issubclass(type(recorded_selection), DoRecorderMixin)


def test_make_recordable():
    # Create a recordable version of NSGA2
    RecordableNSGA2 = make_recordable(NSGA2)

    # Instantiate the recordable NSGA2 with 'selection' to be recorded
    algorithm = RecordableNSGA2(pop_size=92, record=["selection"])

    # Check that the selection operator is now a subclass of DoRecorderMixin
    assert isinstance(algorithm.mating.selection, DoRecorderMixin)
    assert issubclass(type(algorithm.mating.selection), DoRecorderMixin)
