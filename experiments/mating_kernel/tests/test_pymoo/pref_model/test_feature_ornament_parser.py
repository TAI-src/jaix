from mating_kernel.pymoo.pref_model.ornament.feature_parser import (
    FeatureOrnamentParser,
    AttributeType,
    Feature,
)
from pymoo.core.individual import Individual
import numpy as np
import pytest


from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from .. import create_pop

features = [
    Feature(name="X", attribute_type=AttributeType.DISTANCE, xl=0.0),
    Feature(name="Y", attribute_type=AttributeType.DISTANCE, xl=0.0),
    Feature(name="r", attribute_type=AttributeType.VALUE, xl=0.0, xu=5),
]
# TODO: Test for recursive once it works


def test_init():
    parser = FeatureOrnamentParser(features)
    assert parser.n_dim == 3
    assert parser.relative is True
    assert np.all(parser.xl == np.array([0.0, 0.0, 0.0]))
    assert np.all(parser.xu == np.array([np.inf, np.inf, 5.0]))

    parser = FeatureOrnamentParser([features[2]])
    assert parser.n_dim == 1
    assert parser.relative is False


class DummyIndividual(Individual):
    def __init__(self, X, Y, r=None):
        self.X = X
        self.Y = Y
        self.r = r


def test_attr_distance():

    mate = DummyIndividual(X=np.array([1.0, 2.0]), Y=np.array([3.0, 4.0]))
    parents = [
        DummyIndividual(X=np.array([1.0, 1.0]), Y=np.array([3.0, 3.0])),
        DummyIndividual(X=np.array([2.0, 2.0]), Y=np.array([4.0, 4.0])),
    ]

    dist_X = FeatureOrnamentParser.attr_distance("X", mate, parents)
    dist_Y = FeatureOrnamentParser.attr_distance("Y", mate, parents)

    assert np.isclose(dist_X, 1.0)
    assert np.isclose(dist_Y, 1.0)


@pytest.mark.skip(
    reason="This test is not working because the mate does not have the attribute 'data.rank' directly. It is nested inside a dictionary."
)
def test_attr_val():
    mate = Individual()
    mate.data = {"rank": 3}

    val_rank = FeatureOrnamentParser.attr_val("data.rank", mate)
    assert val_rank == 3


def test_ornament():
    mate = DummyIndividual(X=np.array([1.0, 2.0]), Y=np.array([3.0, 4.0]), r=3)
    parents = [
        DummyIndividual(X=np.array([1.0, 1.0]), Y=np.array([3.0, 3.0]), r=2),
        DummyIndividual(X=np.array([2.0, 2.0]), Y=np.array([4.0, 4.0]), r=1),
    ]

    parser = FeatureOrnamentParser(features)
    ornament_values = parser._ornament(mate, parents=parents)

    assert ornament_values.shape == (3,)
    assert np.isclose(ornament_values[0], 1.0)  # Distance for X
    assert np.isclose(ornament_values[1], 1.0)  # Distance for Y
    assert np.isclose(ornament_values[2], 3.0)  # Value for data.rank


def test_pred_feat_ornament_parser():
    problem = REProblem(REProblemConfig(), inst=0)
    features = FeatureOrnamentParser.get_features(
        problem,
        kernel=True,
        abs=False,
        rel_fit=False,
        state=False,
        age=False,
        no_x=False,
    )
    parser = FeatureOrnamentParser(features)
    assert parser.n_dim == 2
    assert parser.relative is True

    pop = create_pop(problem, size=5)
    mate = pop[0]
    parents = [pop[1], pop[2]]
    ornament_values = parser.ornament(
        mate, parents=parents, population=None, problem=None
    )
    assert ornament_values.shape == (2,)
    assert np.all(ornament_values >= parser.xl)
    assert np.all(ornament_values <= parser.xu)
