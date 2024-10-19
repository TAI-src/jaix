import gymnasium as gym
from jaix.env.wrapper import AnyFitWrapper
from gymnasium.utils.env_checker import check_env
from . import DummyEnv


def test_init():
    env = gym.make("MountainCar-v0", render_mode="rgb_array")
    wrapped_env = AnyFitWrapper(env)
    assert wrapped_env.steps == 0


def test_default():
    wrapped_env = AnyFitWrapper(DummyEnv())
    check_env(wrapped_env, skip_render_check=True)


def test_post_log_scale_axis():
    assert AnyFitWrapper._pos_log_scale_axis(10) == 1
    assert AnyFitWrapper._pos_log_scale_axis(100) == 2
    assert AnyFitWrapper._pos_log_scale_axis(0.1) == 0
    assert AnyFitWrapper._pos_log_scale_axis(1) == 0
    assert AnyFitWrapper._pos_log_scale_axis(-1) == 0
    assert AnyFitWrapper._pos_log_scale_axis(-0.5) == 0
    assert AnyFitWrapper._pos_log_scale_axis(-10) == -1


def test_transform_r():
    env = gym.make("MountainCar-v0", render_mode="rgb_array")
    wrapped_env = AnyFitWrapper(env)

    wrapped_env.first_r = 10.5
    wrapped_env.steps = 1
    assert wrapped_env.transform_r(0.5) == 1
    assert wrapped_env.best_r == 0.5
    assert wrapped_env.transform_r(100) == 1
    assert wrapped_env.best_r == 0.5
    assert wrapped_env.transform_r(0.05) >= 1
    assert wrapped_env.best_r == 0.05

    wrapped_env.steps = 10
    wrapped_env.best_r = 5
    assert wrapped_env.transform_r(0.5) == 0.5
    assert wrapped_env.transform_r(0.05) >= 0.5

    wrapped_env.steps = 1
    wrapped_env.first_r = 0.5
    wrapped_env.best_r = 0.5
    assert wrapped_env.transform_r(-9.5) == 1
    wrapped_env.best_r = 0.5
    assert wrapped_env.transform_r(-99.5) == 2


def test_wrappers():
    env = gym.make("MountainCar-v0", render_mode="rgb_array")
    wrapped_env = AnyFitWrapper(env)
    wrapped_env.reset(seed=1337)

    env_cpy = gym.make("MountainCar-v0", render_mode="rgb_array")
    env_cpy.reset(seed=1337)

    num_steps = 1000
    rewards = []
    for i in range(num_steps):
        act = wrapped_env.action_space.sample()
        wrapped_env.step(act)
        assert wrapped_env.steps == i + 1
        _, r, _, _, _ = env_cpy.step(act)
        rewards.append(r)
    assert rewards[0] == wrapped_env.first_r
    assert min(rewards) == wrapped_env.best_r
