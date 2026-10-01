import numpy as np
import pytest
from pymoo.core.individual import Individual
from pymoo.core.population import Population

from mating_kernel.pymoo.parser.population_parser import PopulationParser


@pytest.mark.parametrize("ideal", [None, np.array([0, 0])])
def test_parse(ideal):
    # create pymoo population with individuals having data attributes
    pop = Population.create(
        *[
            Individual(X=[i, 0], F=[i * 2], data={"attr1": i, "attr2": i + 1})
            for i in range(5)
        ]
    )
    simulated_selection_records = [{"args": "a", "kwargs": {}, "output": pop}]

    parser = PopulationParser(ideal=ideal)
    parsed_data = parser.parse({"survival": simulated_selection_records})
    stats = parsed_data[0]
    assert "attr1_mean" in stats
    assert stats["attr1_mean"] == 2.0  # mean of [0, 1, 2, 3, 4]
    assert "attr2_min" in stats
    assert stats["attr2_min"] == 1.0  # min of [1, 2, 3, 4, 5]
    if ideal is not None:
        assert "dist_to_ideal_mean" in stats
    else:
        assert "dist_to_ideal_mean" not in stats
