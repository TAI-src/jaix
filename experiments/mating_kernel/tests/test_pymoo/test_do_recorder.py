import numpy as np
from pymoo.core.population import Population
from pymoo.operators.selection.tournament import TournamentSelection
from pymoo.problems import get_problem

from mating_kernel.pymoo.do_recorder import DoRecorderMixin


class RecordedTournamentSelection(DoRecorderMixin, TournamentSelection):
    pass


def dummy_comp(pop, P, **kwargs):
    return P[:, 0]


def test_do_recorder():
    # Create a TournamentSelection instance
    selection = RecordedTournamentSelection(func_comp=dummy_comp)

    # Create a population of individuals
    problem = get_problem("zdt1")
    x = np.random.uniform(low=problem.xl, high=problem.xu, size=(10, len(problem.xl)))

    pop = Population.new("X", x)
    pop.set("F", problem.evaluate(pop.get("X")))
    off = selection.do(None, pop, n_select=3, n_parents=2)

    assert len(selection.records) == 1
    records = selection.retrieve_records()
    assert np.array_equal(records[0]["output"], off)
    assert len(records) == 1
    assert len(selection.records) == 0

    off2 = selection.do(None, pop, n_select=2, n_parents=2)
    off3 = selection.do(None, pop, n_select=1, n_parents=2)
    assert len(selection.records) == 2
    records = selection.retrieve_records()
    assert np.array_equal(records[0]["output"], off2)
    assert np.array_equal(records[1]["output"], off3)
    assert len(records) == 2
    assert len(selection.records) == 0
