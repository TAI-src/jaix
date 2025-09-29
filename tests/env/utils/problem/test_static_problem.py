from jaix.env.utils.problem.sphere import Sphere, SphereConfig
from jaix.env.utils.problem.static_problem import StaticProblem
from ttex.config import ConfigurableObjectFactory as COF
import pickle
import os


def test_init():
    # Test that defaults are set correctly
    keys = ["lower_bounds", "upper_bounds", "min_values", "max_values"]
    stat_prob = StaticProblem(3, 2)
    for key in keys:
        assert hasattr(stat_prob, key)


def test_evals_left():
    config = SphereConfig(
        dimension=3,
        num_objectives=2,
        mult=1,
        x_shifts=[[0, 0, 0], [0, 0, 0]],
        y_shifts=[0, 0],
        precision=1e-8,
    )
    sphere = COF.create(Sphere, config, 1)
    assert sphere.evalsleft(budget_multiplier=2) == 6
    sphere([0, 0, 0])
    assert sphere.evalsleft(budget_multiplier=2) == 5


def test_stop_evals():
    config = SphereConfig(
        dimension=3,
        num_objectives=2,
        mult=1,
        x_shifts=[[0, 0, 0], [0, 0, 0]],
        y_shifts=[0, 0],
        precision=1e-8,
    )
    sphere = COF.create(Sphere, config, 1)
    sphere([1, 2, 3])
    assert not sphere.stop(1)
    assert not sphere.final_target_hit()
    sphere([1, 2, 3])
    assert not sphere.stop(1)
    assert not sphere.final_target_hit()
    sphere([1, 2, 3])
    assert sphere.stop(1)
    assert not sphere.final_target_hit()


def test_stop_target_hit():
    config = SphereConfig(
        dimension=3,
        num_objectives=2,
        mult=1,
        x_shifts=[[1, 2, 3], [1, 2, 3]],
        y_shifts=[4, 4],
        precision=1e-8,
    )
    sphere = COF.create(Sphere, config, 1)
    sphere([0, 0, 0])
    assert not sphere.stop(1000)
    assert not sphere.final_target_hit()
    sphere([1, 2, 3])
    assert sphere.stop(1000)
    assert sphere.final_target_hit()


"""
def test_recommend_basic():
    config = SphereConfig(
        dimension=3,
        num_objectives=2,
        mult=1,
        x_shifts=[[0, 0, 0], [0, 0, 0]],
        y_shifts=[0, 0],
        precision=1e-8,
    )
    sphere = COF.create(Sphere, config)
    sphere.recommend([0, 0, 0])
    sphere.recommend([1, 1, 1])
    assert sphere.recommendations == {0: [[0, 0, 0], [1, 1, 1]]}
    sphere([3, 3, 3])
    assert sphere.recommendations == {0: [[0, 0, 0], [1, 1, 1]], 1: [[1, 1, 1]]}
    sphere([4, 4, 4])
    assert sphere.recommendations == {0: [[0, 0, 0], [1, 1, 1]], 1: [[1, 1, 1]], 2: []}
    sphere([5, 5, 5])
    assert sphere.recommendations == {
        0: [[0, 0, 0], [1, 1, 1]],
        1: [[1, 1, 1]],
        2: [[5, 5, 5]],
    }
    sphere.recommend([2, 2, 2])
    assert sphere.recommendations == {
        0: [[0, 0, 0], [1, 1, 1]],
        1: [[1, 1, 1]],
        2: [[5, 5, 5]],
        3: [[2, 2, 2]],
    }
    sphere.recommend([-1, -1, -1])
    sphere.recommend([4, 4, 4])
    assert sphere.recommendations == {
        0: [[0, 0, 0], [1, 1, 1]],
        1: [[1, 1, 1]],
        2: [[5, 5, 5]],
        3: [[2, 2, 2], [-1, -1, -1], [4, 4, 4]],
    }
    sphere([-2, -2, -2])
    sphere.recommend([-3, -3, -3])
    assert sphere.recommendations == {
        0: [[0, 0, 0], [1, 1, 1]],
        1: [[1, 1, 1]],
        2: [[5, 5, 5]],
        3: [[2, 2, 2], [-1, -1, -1], [4, 4, 4]],
        4: [[-3, -3, -3]],
    }
"""


def test_close():
    config = SphereConfig(
        dimension=3,
        num_objectives=2,
        mult=1,
        x_shifts=[[0, 0, 0], [0, 0, 0]],
        y_shifts=[0, 0],
        precision=1e-8,
    )
    sphere = COF.create(Sphere, config, 1)
    sphere.recommend([0, 0, 0])
    rec_file_name = sphere.close()
    with open(rec_file_name, "rb") as rec_file:
        recs = pickle.load(rec_file)
    assert recs == {0: [[0, 0, 0]]}
    # clean up file
    os.remove(rec_file_name)
