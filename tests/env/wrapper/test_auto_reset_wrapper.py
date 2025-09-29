from jaix.env.wrapper.auto_reset_wrapper import AutoResetWrapper, AutoResetWrapperConfig
import gymnasium as gym
from . import DummyEnv
import pytest
from gymnasium.utils.env_checker import check_env


def test_init():
    env = DummyEnv()
    config = AutoResetWrapperConfig(min_steps=2)
    wrapped_env = AutoResetWrapper(config, env)
    assert wrapped_env.steps == 0
    assert wrapped_env.failed_resets == 0
    assert wrapped_env.min_steps == 2


def test_default():
    wrapped_env = AutoResetWrapper(AutoResetWrapperConfig(min_steps=2), DummyEnv())
    check_env(wrapped_env, skip_render_check=True)


def test_normal_exec():
    env = DummyEnv()
    config = AutoResetWrapperConfig(min_steps=2)
    wrapped_env = AutoResetWrapper(config, env)
    _, info = wrapped_env.reset()
    assert "final_r" not in info
    assert wrapped_env.man_resets == 1
    assert wrapped_env.steps == 0
    assert wrapped_env.prev_r == None

    _, r, term, trunc, info = wrapped_env.step(wrapped_env.action_space.sample())
    assert not trunc
    assert not term
    assert wrapped_env.prev_r == r
    assert "final_info" not in info
    assert "final_r" not in info

    assert not wrapped_env.stop()


@pytest.mark.parametrize("num_steps", [1, 2, 3])
def test_final_r_step(num_steps):
    min_steps = 3
    env = DummyEnv()
    config = AutoResetWrapperConfig(min_steps=min_steps)
    wrapped_env = AutoResetWrapper(config, env)
    _, info = wrapped_env.reset()

    r = None
    for _ in range(num_steps):
        _, r, _, _, info = wrapped_env.step(wrapped_env.action_space.sample())
        assert "final_info" not in info

    wrapped_env.unwrapped._term = True
    _, _, _, _, info = wrapped_env.step(wrapped_env.action_space.sample())

    assert "final_info" in info
    if num_steps >= min_steps:
        assert info["final_r"] == r
    else:
        assert "final_r" not in info
    assert wrapped_env.auto_resets == 1
    assert wrapped_env.steps == 1
    # Make sure that the inner environment got reset
    assert not wrapped_env.unwrapped._term

    if num_steps >= min_steps:
        assert not wrapped_env.stop()
    else:
        assert wrapped_env.stop()


@pytest.mark.parametrize("num_steps", [1, 2, 3])
def test_final_r_reset(num_steps):
    min_steps = 3
    env = DummyEnv()
    config = AutoResetWrapperConfig(min_steps=min_steps)
    wrapped_env = AutoResetWrapper(config, env)
    _, info = wrapped_env.reset()

    prev_r = None
    for _ in range(num_steps):
        _, prev_r, _, _, info = wrapped_env.step(wrapped_env.action_space.sample())
        assert "final_info" not in info

    wrapped_env.unwrapped._trunc = True
    _, info = wrapped_env.reset()
    assert info["final_r"] == prev_r
    assert wrapped_env.man_resets == 2
    assert wrapped_env.steps == 0
    assert not wrapped_env.unwrapped._trunc

    assert not wrapped_env.stop()


@pytest.mark.parametrize("failed_resets", [1, 2, 3])
def test_stop(failed_resets):
    env = DummyEnv()
    config = AutoResetWrapperConfig(failed_resets_thresh=failed_resets)
    wrapped_env = AutoResetWrapper(config, env)

    _, _ = wrapped_env.reset()

    for i in range(1, failed_resets + 1):
        assert not wrapped_env.stop()
        wrapped_env.reset()
        wrapped_env.unwrapped._term = True
        wrapped_env.step(wrapped_env.action_space.sample())
        assert wrapped_env.failed_resets == i
        assert wrapped_env.steps == 1
        assert wrapped_env.auto_resets == i
        assert not wrapped_env.unwrapped._term

    assert wrapped_env.stop()
