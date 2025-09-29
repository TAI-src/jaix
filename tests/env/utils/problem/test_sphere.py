from jaix.env.utils.problem.sphere import Sphere, SphereConfig
from ttex.config import ConfigurableObjectFactory as COF
import pytest


def test_init():
    config = SphereConfig(
        dimension=3,
        num_objectives=2,
        mult=1,
        x_shifts=[[0, 0, 0], [0, 0, 0]],
        y_shifts=[0, 0],
        precision=1e-8,
    )
    sphere = COF.create(Sphere, config, 1)
    for key in ["dimension", "mult", "precision", "lower_bounds", "upper_bounds"]:
        assert getattr(sphere, key) == getattr(config, key)
    assert not sphere.final_target_hit()


def test_final_target_hit():
    config = SphereConfig(
        dimension=3,
        num_objectives=2,
        mult=1,
        x_shifts=[[0, 0, 0], [0, 0, 0]],
        y_shifts=[0, 0],
        precision=1e-8,
    )
    sphere = COF.create(Sphere, config, 1)
    f = sphere([0, 0, 1e-15])  # Slightly offset to check precision works
    assert sphere.final_target_hit()
    config.x_shifts = [[0, 0, 0], [1, 2, 3]]
    sphere = COF.create(Sphere, config, 1)
    assert not sphere.final_target_hit()
    f = sphere([0, 0, 1e-15])  # Slightly offset to check precision works
    assert not sphere.final_target_hit()


def test_eval():
    config = SphereConfig(
        dimension=3,
        num_objectives=2,
        mult=1,
        x_shifts=[[1, 2, 3], [1, 2, 3]],
        y_shifts=[4, 4],
        precision=1e-8,
    )
    sphere = COF.create(Sphere, config, 1)
    f, _ = sphere([1, 2, 3])
    assert sphere.final_target_hit()
    assert f == [4, 4]


@pytest.mark.parametrize("n_obj", [1, 2])
def test_format(n_obj):
    config = SphereConfig(
        dimension=3,
        num_objectives=n_obj,
        mult=1,
        x_shifts=[[1, 2, 3]] * n_obj,
        y_shifts=[4] * n_obj,
        precision=1e-8,
    )
    sphere = COF.create(Sphere, config, 1)
    f, _ = sphere([1, 2, 3])
    assert len(f) == n_obj
