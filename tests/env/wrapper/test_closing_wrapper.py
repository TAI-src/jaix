from jaix.env.wrapper.closing_wrapper import ClosingWrapper
from . import DummyEnv
from gymnasium.utils.env_checker import check_env
import pytest


def test_init():
    env = DummyEnv()
    kwargs = {}
    wrapped_env = ClosingWrapper(env, **kwargs)
    assert not wrapped_env.closed


def test_default():
    wrapped_env = ClosingWrapper(DummyEnv())
    check_env(wrapped_env, skip_render_check=True)


def test_closed():
    wrapped_env = ClosingWrapper(DummyEnv())
    wrapped_env.close()
    assert wrapped_env.closed

    with pytest.raises(ValueError):
        wrapped_env.reset(seed=1)
    with pytest.raises(ValueError):
        wrapped_env.step(wrapped_env.action_space.sample())
