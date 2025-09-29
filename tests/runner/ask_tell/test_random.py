from jaix.runner.ask_tell.strategy.random import RandomATStrat, RandomATStratConfig
from . import DummyEnv, loop
import pytest


def test_init():
    env = DummyEnv()
    rand_opt = RandomATStrat(
        RandomATStratConfig(ask_size=5), [env.action_space.sample()]
    )
    assert env.action_space.contains(rand_opt.xcurrent[0])


@pytest.mark.parametrize("dimension,num_objectives", [(3, 1), (3, 2)])
def test_loop(dimension, num_objectives):
    rand_opt = RandomATStrat(RandomATStratConfig(ask_size=5), [[0] * dimension])

    X, Y = loop(dimension, num_objectives, rand_opt)
    assert len(X) == 5  # ask size
