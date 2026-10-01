from mating_kernel.pymoo.parser.reproduction_parser import ReproductionParser
from pymoo.core.individual import Individual
from pymoo.operators.crossover.sbx import SBX
import numpy as np
import pytest


class DummyProblem:
    def __init__(self):
        self.n_var = 2
        self.xl = np.array([0, 0])
        self.xu = np.array([200, 200])


@pytest.mark.parametrize("n_offsprings", [2, 3])
def test_extract_lineage(n_offsprings):
    if n_offsprings != 2:
        pytest.skip("SBX crossover only supports 2 offspring")
    # create parent pairings in a way that it will be obvious which offspring came from which parents
    parent1 = Individual(X=[0, 0], F=[0], data={"id": "p1"})
    parent2 = Individual(X=[1, 1], F=[1], data={"id": "p2"})
    parent3 = Individual(X=[100, 100], F=[2], data={"id": "p3"})
    parent4 = Individual(X=[101, 101], F=[3], data={"id": "p4"})
    parents = [[parent1, parent2], [parent3, parent4]]  # two mating events
    # actually do a crossover to generate offspring
    crossover = SBX(prob=1.0, eta=15, n_offsprings=n_offsprings)
    offspring = crossover.do(
        DummyProblem(), parents
    )  # two or 3 offspring per mating events
    lineage = ReproductionParser.extract_lineage(
        parents, offspring, n_offsprings=n_offsprings
    )
    assert (
        len(lineage) == len(parents) * n_offsprings
    )  # 2 mating events * 2 offspring per event
    # check that the parents in the lineage match the original parents
    for i, entry in enumerate(lineage):
        offspring = entry["offspring"].X
        parents = [p.X for p in entry["parents"]]

        for j in range(len(offspring)):
            min_val = min(parents[0][j], parents[1][j])
            max_val = max(parents[0][j], parents[1][j])
            tolerance = 10  # SBX can produce values outside the range of parents, but should be within a reasonable range
            assert min_val - tolerance <= offspring[j] <= max_val + tolerance


def test_parse_individual():
    ind = Individual(X=[1, 2], F=[3], data={"attr": 4})
    parser = ReproductionParser()
    parsed = parser.parse_individual(ind)
    assert parsed["X"] == [1, 2]
    assert parsed["F"] == [3]
    assert parsed["attr"] == 4


def test_parse():
    # create parent pairings in a way that it will be obvious which offspring came from which parents
    parent1 = Individual(X=[0, 0], F=[0], data={"id": "p1"})
    parent2 = Individual(X=[1, 1], F=[1], data={"id": "p2"})
    parent3 = Individual(X=[100, 100], F=[2], data={"id": "p3"})
    parent4 = Individual(X=[101, 101], F=[3], data={"id": "p4"})
    parents = [[parent1, parent2], [parent3, parent4]]  # two mating events
    # actually do a crossover to generate offspring
    crossover = SBX(prob=1.0, eta=15, n_offsprings=2)
    offspring = crossover.do(DummyProblem(), parents)  # two offspring per mating events
    survived = [
        offspring[0],
        offspring[-1],
    ]  # only the first and last offspring survived

    # simulate the data structure that would be passed to the parser
    data = {
        "mating.selection": [{}, {"args": None, "kwargs": None, "output": parents}],
        "mating.crossover": [
            {"n_offsprings": 2},
            {"args": None, "kwargs": None, "output": offspring},
        ],
        "survival": [{}, {"args": None, "kwargs": None, "output": survived}],
    }

    parser = ReproductionParser()
    parsed_data = parser.parse(data)

    assert len(parsed_data) == len(offspring)
    for entry in parsed_data:
        assert "o_X" in entry
        assert "o_F" in entry
        assert "p0_X" in entry
        assert "p0_F" in entry
        assert "p1_X" in entry
        assert "p1_F" in entry
        assert "survived" in entry
    assert parsed_data[0]["survived"] is True
    assert parsed_data[1]["survived"] is False
    assert parsed_data[2]["survived"] is False
    assert parsed_data[3]["survived"] is True
