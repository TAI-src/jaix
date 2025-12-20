import gymnasium as gym
from jaix.env.wrapper.improvement_reward_wrapper import (
    ImprovementRewardWrapper,
    ImprovementRewardWrapperConfig,
    ImprovementType,
)
from gymnasium.utils.env_checker import check_env
from . import DummyEnv
import pytest


def test_init():
    env = gym.make("MountainCar-v0", render_mode="rgb_array")
    config = ImprovementRewardWrapperConfig()
    wrapped_env = ImprovementRewardWrapper(config, env)
    assert wrapped_env.imp_type == config.imp_type
    assert wrapped_env.state_eval == config.state_eval


def test_default():
    config = ImprovementRewardWrapperConfig()
    wrapped_env = ImprovementRewardWrapper(config, DummyEnv())
    check_env(wrapped_env, skip_render_check=True)


def test_post_log_scale_axis():
    assert ImprovementRewardWrapper._pos_log_scale_axis(10) == 1
    assert ImprovementRewardWrapper._pos_log_scale_axis(100) == 2
    assert ImprovementRewardWrapper._pos_log_scale_axis(0.1) == 0
    assert ImprovementRewardWrapper._pos_log_scale_axis(1) == 0
    assert ImprovementRewardWrapper._pos_log_scale_axis(-1) == 0
    assert ImprovementRewardWrapper._pos_log_scale_axis(-0.5) == 0
    assert ImprovementRewardWrapper._pos_log_scale_axis(-10) == -1


@pytest.mark.parametrize("min_bool", [True, False])
def test_compute_imp(min_bool):
    env = gym.make("MountainCar-v0", render_mode="rgb_array")
    config = ImprovementRewardWrapperConfig(is_min=min_bool, transform=True)
    wrapped_env = ImprovementRewardWrapper(config, env)

    assert wrapped_env._compute_imp(10, 5) == (5 if min_bool else 0)
    assert wrapped_env._compute_imp(5, 10) == (0 if min_bool else 5)
    assert wrapped_env._compute_imp(5, 5) == 0


def test__get_improvement():
    env = gym.make("MountainCar-v0", render_mode="rgb_array")
    config = ImprovementRewardWrapperConfig(
        transform=False, imp_type=ImprovementType.OVER_FIRST, is_min=True
    )
    wrapped_env = ImprovementRewardWrapper(config, env)
    assert wrapped_env.first_val is None
    assert wrapped_env.best_val is None
    assert wrapped_env.last_val is None
    wrapped_env.steps = 1
    assert wrapped_env._get_improvement(0) == 0.0

    wrapped_env.first_val = 10.5
    wrapped_env.best_val = 5
    wrapped_env.last_val = 7
    wrapped_env.steps = 2
    assert wrapped_env._get_improvement(0.5) == 10.0
    assert wrapped_env._get_improvement(11) == 0

    wrapped_env.imp_type = ImprovementType.OVER_BEST
    assert wrapped_env._get_improvement(0.5) == 4.5
    assert wrapped_env._get_improvement(6) == 0

    wrapped_env.imp_type = ImprovementType.OVER_LAST
    assert wrapped_env._get_improvement(0.5) == 6.5
    assert wrapped_env._get_improvement(8) == 0

    wrapped_env.imp_type = ImprovementType.BEST_SINCE_FIRST
    assert wrapped_env._get_improvement(0.5) == 10
    assert wrapped_env._get_improvement(7) == 5.5


def test_get_improvement():
    env = gym.make("MountainCar-v0", render_mode="rgb_array")
    config = ImprovementRewardWrapperConfig(
        transform=True,
        imp_type=ImprovementType.BEST_SINCE_FIRST,
        state_eval="r",
        is_min=True,
    )
    wrapped_env = ImprovementRewardWrapper(config, env)

    wrapped_env.best_val = 100000
    wrapped_env.first_val = 10.5
    wrapped_env.last_val = 5000
    wrapped_env.steps = 0
    assert wrapped_env.get_improvement(0.5) == 1

    wrapped_env.best_val = 0.5
    assert wrapped_env.get_improvement(100) == 1
    assert wrapped_env.get_improvement(0.05) >= 1

    wrapped_env.steps = 9
    wrapped_env.best_val = 5
    assert wrapped_env.get_improvement(0.5) == 0.5
    assert wrapped_env.get_improvement(0.05) >= 0.5

    wrapped_env.steps = 0
    wrapped_env.first_val = 0.5
    wrapped_env.best_val = 0.5
    assert wrapped_env.get_improvement(-9.5) == 1
    assert wrapped_env.get_improvement(-99.5) == 2


@pytest.mark.parametrize("state_eval", ["obs0", "r", "val"])
def test_val_choice(state_eval):
    env = DummyEnv(num_objectives=1)
    config = ImprovementRewardWrapperConfig(state_eval=state_eval)
    wrapped_env = ImprovementRewardWrapper(config, env)
    wrapped_env.reset()
    obs, r, _, _, info = wrapped_env.step(wrapped_env.action_space.sample())
    obs, r, _, _, info = wrapped_env.step(wrapped_env.action_space.sample())
    if state_eval == "obs0":
        assert wrapped_env.last_val == obs[0]
    elif state_eval == "r":
        assert wrapped_env.last_val != r
    elif state_eval == "val":
        assert wrapped_env.last_val == info["val"]

    assert info[f"raw_{state_eval}"] == wrapped_env.last_val
