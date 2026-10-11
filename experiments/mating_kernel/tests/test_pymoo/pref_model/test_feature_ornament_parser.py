from mating_kernel.pymoo.pref_model.ornament.feature_parser import (
    FeatureOrnamentParser,
    AttributeType,
    Feature,
    ind_getattr,
    recursive_getdict,
)
from pymoo.core.individual import Individual
import numpy as np
import pytest


from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from .. import create_pop


def test_recursive_getdict():
    d = {"a": {"b": {"c": 1}}}
    assert recursive_getdict(d, ["a", "b", "c"]) == 1
    assert recursive_getdict(d, ["a", "b"]) == {"c": 1}
    assert recursive_getdict(d, ["a", "x"], default="not found") == "not found"
    assert recursive_getdict(d, ["a", "b", "c", "d"], default=None) is None


def test_ind_getattr():
    ind = Individual()
    ind.X = np.array([1.0, 2.0])
    ind.data = {"rank": 3, "nested": {"value": 5}}

    assert np.array_equal(ind_getattr(ind, "X"), np.array([1.0, 2.0]))
    assert ind_getattr(ind, "data.rank") == 3
    assert ind_getattr(ind, "data.nested.value") == 5
    assert np.isnan(ind_getattr(ind, "data.nonexistent"))


features = [
    Feature(
        name="X",
        attribute_type=AttributeType.DISTANCE,
        xl=np.array([0]),
        xu=np.array([np.inf]),
    ),
    Feature(
        name="Y",
        attribute_type=AttributeType.DISTANCE,
        xl=np.array([0]),
        xu=np.array([np.inf]),
    ),
    Feature(
        name="data.rank",
        attribute_type=AttributeType.VALUE,
        xl=np.array([0]),
        xu=np.array([5]),
    ),
    Feature(
        name="Y",
        attribute_type=AttributeType.VALUE,
        xl=np.array([-np.inf, -np.inf]),
        xu=np.array([np.inf, np.inf]),
    ),
]


def test_init():
    parser = FeatureOrnamentParser(features)
    assert parser.n_dim == 5  # 2 distances + 1 value + 2 values
    assert parser.relative is True
    assert np.all(parser.xl == np.array([0.0, 0.0, 0.0, -np.inf, -np.inf]))
    assert np.all(parser.xu == np.array([np.inf, np.inf, 5.0, np.inf, np.inf]))

    parser = FeatureOrnamentParser([features[2]])
    assert parser.n_dim == 1
    assert parser.relative is False


class DummyIndividual(Individual):
    def __init__(self, X, Y, r=None):
        self.X = X
        self.Y = Y
        self.data = {"rank": r} if r is not None else {}


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

    assert ornament_values.shape == (5,)
    assert np.isclose(ornament_values[0], 1.0)  # Distance for X
    assert np.isclose(ornament_values[1], 1.0)  # Distance for Y
    assert np.isclose(ornament_values[2], 3.0)  # Value for data.rank


@pytest.mark.parametrize("no_x", [True, False])
def test_get_features(no_x):
    problem = REProblem(REProblemConfig(), inst=0)
    features = FeatureOrnamentParser.get_features(
        problem,
        kernel=True,
        abs=True,
        rel_fit=True,
        state=False,
        age=True,
        no_x=no_x,
    )
    expected_length = 2 + problem.dimension + problem.num_objectives + 2 + 1
    if no_x:
        expected_length -= 1 + problem.dimension
    expected_num_features = 2 + 2 + 2 + 1
    if no_x:
        expected_num_features -= 1 + 1
    assert len(features) == expected_num_features

    parser = FeatureOrnamentParser(features)
    assert parser.n_dim == expected_length


@pytest.mark.parametrize("nan", [True, False])
def test_pred_feat_ornament_parser(nan):
    problem = REProblem(REProblemConfig(), inst=0)
    features = FeatureOrnamentParser.get_features(
        problem,
        kernel=True,
        abs=True,
        rel_fit=True,
        state=False,
        age=True,
        no_x=False,
    )
    parser = FeatureOrnamentParser(features)
    assert parser.n_dim == 2 + problem.dimension + problem.num_objectives + 2 + 1
    assert parser.relative is True

    pop = create_pop(problem, size=5)
    for ind in pop:
        ind.data["n_gen"] = np.random.randint(0, 10)
        ind.data["crowding"] = np.random.uniform(0, 1)
        ind.data["rank"] = np.random.randint(0, 5)
    if nan:
        pop[0].data["crowding"] = np.nan  # Introduce a NaN value for testing
        pop[2].data["n_gen"] = np.nan  # Introduce a NaN value for testing
    mate = pop[0]
    parents = [pop[1], pop[2]]
    ornament_values = parser.ornament(
        mate, parents=parents, population=None, problem=None
    )
    assert ornament_values.shape == (parser.n_dim,)
    for i, (o_val, xl, xu) in enumerate(zip(ornament_values, parser.xl, parser.xu)):
        if (
            i
            in [
                9,
                10,
            ]
            and nan
        ):  # These are the nan values we introduced, so we skip the assertion for them
            continue
        assert not np.isnan(o_val), f"Value at {i} is NaN"
        if i not in [
            2,
            3,
        ]:  # Skip F boundary tests (since we are using nadir, not max
            assert (
                xl <= o_val <= xu
            ), f"Value {o_val} at {i} is out of bounds [{xl}, {xu}]"
