from jaix.runner.ask_tell import ATOptimiserConfig, ATOptimiser
from jaix.runner.ask_tell.strategy import RandomATStrat, RandomATStratConfig
from ttex.config import ConfigurableObjectFactory as COF
from . import DummyEnv, loop
import pytest


def create_at_opt(dimension, num_objectives):
    env = DummyEnv(dimension, num_objectives)
    config = ATOptimiserConfig(
        strategy_class=RandomATStrat,
        strategy_config=RandomATStratConfig(ask_size=5),
        init_pop_size=1,
        stop_after=3,
    )
    at_opt = COF.create(ATOptimiser, config, env)
    return at_opt


def test_init():
    at_opt = create_at_opt(6, 2)
    assert at_opt.countiter == 0
    assert at_opt.name == "Random Search"


@pytest.mark.parametrize("dimension,num_objectives", [(3, 1), (3, 2)])
def test_loop(dimension, num_objectives):
    at_opt = create_at_opt(dimension, num_objectives)

    X, Y = loop(dimension, num_objectives, at_opt)
    assert len(X) == 5  # ask size


def test_stop():
    at_opt = create_at_opt(6, 2)
    for _ in range(at_opt.stop_after):
        X, Y = loop(6, 2, at_opt)
    assert at_opt.stop()["countiter"] == at_opt.stop_after