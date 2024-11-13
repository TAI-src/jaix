from jaix.env.utils.problem.rbf import RBFAdapterConfig, RBFAdapter
import pytest
from math import isclose
import numpy as np


@pytest.mark.parametrize(
    "start,length,num_splits", [(0, 10, 5), (-3, 7.5, 8), (1, 2, 2)]
)
def test_split_range(start, length, num_splits):
    points = RBFAdapter._split_range(start, length, num_splits)
    assert len(points) == num_splits
    assert points[-1] - points[0] == length
    assert points[0] == start

    d = points[1] - points[0]
    assert all([isclose(points[i + 1] - points[i], d) for i in range(num_splits - 1)])


@pytest.mark.parametrize("start,length,num_splits", [(-5, 10, 1)])
def test_split_range_edge(start, length, num_splits):
    points = RBFAdapter._split_range(start, length, num_splits)
    assert len(points) == num_splits


def get_config():
    config = RBFAdapterConfig(
        num_rad=5,
        const_ratio_x=0.5,
        num_measure_points=20,
    )
    return config


@pytest.mark.parametrize("seed", [42, 1337])
def test_setup(seed):
    np.random.seed(seed)
    config = get_config()
    targets, centers = RBFAdapter._setup(config)

    x_length = config.x_val_range[1] - config.x_val_range[0]
    assert len(centers) == config.num_rad
    assert centers[-1] - centers[0] == x_length

    assert len(targets) == config.num_measure_points
    p = [(m, v) for (m, v) in targets if v != 0]
    const_ratio_x = (p[-1][0] - p[0][0]) / x_length
    abs(const_ratio_x - config.const_ratio_x) <= 1 / config.num_measure_points
    assert p[0][1] >= config.y_val_range[0]
    assert p[0][1] <= config.y_val_range[1]


def test_init():
    rbf_adapter1 = RBFAdapter(get_config(), 5)
    rbf_adapter2 = RBFAdapter(get_config(), 5)
    assert rbf_adapter1.targets == rbf_adapter2.targets
    assert rbf_adapter1.centers == rbf_adapter2.centers


def test_comp_fit():
    config = RBFAdapterConfig(
        num_rad=1,
        const_ratio_x=1,
        num_measure_points=1,
    )
    rbf_adapter = RBFAdapter(config, 3)
    # Just check that I can execute and is correct format
    fit = rbf_adapter.comp_fit([0] * (2 * rbf_adapter.num_rad))
    assert isinstance(fit, float)

    # Check value makes sense
    rbf_adapter.err = lambda d: d
    d = rbf_adapter.comp_fit([1, 1])
    rbf = d[0] + rbf_adapter.targets[0][1]
    assert rbf == 1
