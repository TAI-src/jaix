from jaix.env.wrapper.max_eval_wrapper import MaxEvalWrapper, MaxEvalWrapperConfig
from . import DummyWrapper, DummyEnv
import pytest


@pytest.fixture(scope="function")
def wenv():
    env = DummyEnv()
    config = MaxEvalWrapperConfig(max_evals=10)
    wenv = MaxEvalWrapper(config, env)
    return wenv


def test_step_stop(wenv):
    assert wenv._evals == 0
    counter = 0

    while not wenv._stop():
        wenv.step(wenv.unwrapped.action_space.sample())
        counter += 1
        assert wenv._evals == counter

    assert counter == 10
    assert wenv._stop()["max_evals"] == 10


def test_reset(wenv):
    wenv._evals = 10
    assert wenv._stop()

    # online
    wenv.reset(options={"online": True})
    assert wenv._stop()

    # not online
    wenv.reset()
    assert not wenv._stop()
    assert wenv._evals == 0
