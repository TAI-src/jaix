from jaix.env.utils.problem.rbf.rbf_adapter import RBFAdapterConfig, RBFAdapter
from jaix.env.utils.problem.rbf.rbf import RBF
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
        num_rad_range=(1, 20),
        ratio_x_range=(0.5, 0.5),
    )
    return config


@pytest.mark.parametrize("seed", [42, 1337])
def test_setup(seed):
    config = get_config()
    rbf_adapter = RBFAdapter(config, seed)

    x_length = config.x_val_range[1] - config.x_val_range[0]
    assert len(rbf_adapter.centers) in range(
        config.num_rad_range[0], config.num_rad_range[1]
    )
    assert rbf_adapter.centers[-1] - rbf_adapter.centers[0] == x_length
    assert rbf_adapter.target_val <= config.y_val_range[1]
    assert rbf_adapter.target_val >= config.y_val_range[0]


@pytest.mark.parametrize("noisy", [True, False])
def test_get_targets(noisy):
    config = get_config()
    config.noisy = noisy
    rbf_adapter = RBFAdapter(config, 1337)

    targets = rbf_adapter.get_targets(config.num_measure_points)

    assert len(targets) == config.num_measure_points
    p = [(m, v) for (m, v) in targets if v != 0]
    # Check that first value point (measure point with non-zero value) is around box start
    x_length = config.x_val_range[1] - config.x_val_range[0]
    tol = x_length / config.num_measure_points * 5
    assert abs(rbf_adapter.box_start - p[0][0]) <= tol
    # Same for last value point
    assert abs(rbf_adapter.box_end - p[-1][0]) <= tol
    # Check that the target value is assigned correctly
    assert p[0][1] == rbf_adapter.target_val

    # _Test noise
    targets_2 = rbf_adapter.get_targets(config.num_measure_points)
    if not noisy:
        assert targets == targets_2
    else:
        assert not targets == targets_2


def test_init():
    rbf_adapter1 = RBFAdapter(get_config(), 5)
    rbf_adapter2 = RBFAdapter(get_config(), 5)
    assert rbf_adapter1.centers == rbf_adapter2.centers


def test_noise():
    config = get_config()
    rbf_adapter = RBFAdapter(config, 1337)
    config.num_measure_points = 10
    fit, r = rbf_adapter.comp_fit([0] * rbf_adapter.num_rad)
    assert not np.isclose(fit, r)


def test_comp_fit():
    config = RBFAdapterConfig(
        num_rad_range=(1, 1),
        ratio_x_range=(1, 1),
        num_measure_points=1,
        noisy=False,
    )
    rbf_adapter = RBFAdapter(config, 3)
    # Just check that I can execute and is correct format
    fit, r = rbf_adapter.comp_fit([0] * rbf_adapter.num_rad)
    assert isinstance(fit, float)
    assert np.isclose(r, fit)

    # Check value makes sense
    rbf_adapter.err = lambda d: d
    d, r = rbf_adapter.comp_fit([1])
    targets = rbf_adapter.get_targets(config.num_measure_points)
    rbf = RBF(rbf_adapter.centers, [1], [1], rbf_adapter.kernel)
    assert d[0] == rbf.eval(targets[0][0]) - targets[0][1]
