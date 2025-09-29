from ttex.config import (
    ConfigurableObjectFactory as COF,
)
from jaix.env.composite.switching_environment import (
    SwitchingEnvironmentConfig,
    SwitchingEnvironment,
)
from jaix.env.utils.switching_pattern.switching_pattern import (
    SeqRegSwitchingPatternConfig,
    SeqRegSwitchingPattern,
)
import pytest
import gymnasium as gym
from gymnasium import spaces
import numpy as np
from jaix.env.utils.problem.sphere import Sphere, SphereConfig
from jaix.env.utils.problem.rbf_fit import RBFFit, RBFFitConfig
from ..utils.problem.rbf.test_rbf_adapter import get_config
from jaix.env.singular.ec_env import (
    ECEnvironment,
    ECEnvironmentConfig,
)
import os
from jaix.env.wrapper.closing_wrapper import ClosingWrapper


@pytest.fixture(scope="function")
def env():
    sp_config = SeqRegSwitchingPatternConfig(wait_period=3)

    env_list = [
        gym.make("MountainCar-v0", render_mode="rgb_array"),
        gym.make("MountainCar-v0", render_mode="rgb_array"),
        gym.make("MountainCar-v0", render_mode="rgb_array"),
    ]

    config = SwitchingEnvironmentConfig(
        SeqRegSwitchingPattern, sp_config, real_time=False
    )
    env = COF.create(
        SwitchingEnvironment,
        config,
        env_list,
    )
    env = ClosingWrapper(env)

    yield env
    files = env.close()
    [os.remove(rec_file) for rec_file in files if rec_file is not None]
    assert env.closed


def test_init(env):
    assert env._current_env == 0
    assert len(env.env_list) == 3
    assert env.steps_counter == 0

    act = env.action_space.sample()
    assert env.env_list[0].action_space.contains(act)

    obs = env.observation_space.sample()
    assert spaces.Discrete(3).contains(obs[0])
    assert env.env_list[0].observation_space.contains(obs[1])


def test_reset(env):
    obs, info = env.reset()

    assert env._current_env == 0
    assert env._timer == 0

    assert obs[0] == 0

    assert "steps_counter" in info["meta"]
    assert all(obs[1] == np.array(env.env_list[0].unwrapped.state, dtype=np.float32))


def intended_change(old, env, idx):
    new = [c_env.unwrapped.state for c_env in env.env_list]
    res = [all(np.array(o) == np.array(n)) for o, n in zip(old, new)]
    assert sum(res) == len(res) - 1  # only one changed
    assert not res[idx]  # correct one changed
    return new


def test_step(env):
    obs, info = env.reset()
    act = env.action_space.sample()
    base_line_states = [c_env.unwrapped.state for c_env in env.env_list]

    obs, _, _, _, info = env.step(act)

    assert obs[0] == 0
    assert info["meta"]["steps_counter"] == 1
    assert info["meta"]["timer"] == 1
    base_line_states = intended_change(base_line_states, env, 0)

    obs, _, _, _, _ = env.step(act)
    assert obs[0] == 0
    base_line_states = intended_change(base_line_states, env, 0)

    obs, _, _, _, _ = env.step(act)
    assert obs[0] == 0
    base_line_states = intended_change(base_line_states, env, 0)

    obs, _, _, _, _ = env.step(act)
    assert obs[0] == 1
    base_line_states = intended_change(base_line_states, env, 1)

    env.reset()
    obs, _, _, _, info = env.step(act)
    assert obs[0] == 0
    assert info["meta"]["steps_counter"] == 5
    assert info["meta"]["timer"] == 1


def test_truncate(env):
    env.steps_counter = 0
    obs, info = env.reset()
    act = env.action_space.sample()

    trunc = False
    while not trunc:
        _, _, _, trunc, info = env.step(act)
    assert info["meta"]["steps_counter"] == 10


def test_render(env):
    env.reset()
    with pytest.raises(gym.error.DependencyNotInstalled):
        env.render()


def test_passthrough(env):
    env.reset()
    # Make sure no errors occurr
    env.steps_counter
    env.unwrapped.state
    env.unwrapped.get_keys_to_action()

    with pytest.raises(AttributeError):
        env.not_exist()


def test_decorator_update(env):
    env.reset()
    env.unwrapped._timer = 100
    act = env.action_space.sample()
    _, _, _, trunc, _ = env.step(act)
    assert trunc


@pytest.fixture(scope="function", params=[Sphere, RBFFit])
def ec_env(request):
    if request.param == Sphere:
        func_config = SphereConfig(
            dimension=3,
            num_objectives=2,
            mult=1,
            x_shifts=[[0, 0, 0], [0, 0, 0]],
            y_shifts=[1, 0],
            precision=1e-8,
        )
    elif request.param == RBFFit:
        func_config = RBFFitConfig(get_config(), 1e-8)
    else:
        raise ValueError()

    n_envs = 3
    env_list = [None] * n_envs
    for i in range(n_envs):
        func = COF.create(request.param, func_config, i)

        config = ECEnvironmentConfig(budget_multiplier=1)
        env = COF.create(ECEnvironment, config, func)
        env_list[i] = env

    sp_config = SeqRegSwitchingPatternConfig(wait_period=5)

    config = SwitchingEnvironmentConfig(
        SeqRegSwitchingPattern, sp_config, real_time=False
    )
    env = COF.create(
        SwitchingEnvironment,
        config,
        env_list,
    )
    env = ClosingWrapper(env)
    yield env
    files = env.close()
    [os.remove(rec_file) for rec_file in files if rec_file is not None]
    assert env.closed


def test_step_through_stop(ec_env):
    ec_env.reset(options={"online": True})
    while not ec_env.stop():
        act = ec_env.action_space.sample()
        ec_env.step(act)


def test_force_stop_time(ec_env):
    if not isinstance(ec_env.env_list[0].func, Sphere):
        pytest.skip("only sphere due to known optimum")
    ec_env.reset(options={"online": True})
    while not ec_env.stop():
        act = [0, 0, 1]  # This will not stop the env
        _, _, _, _, info = ec_env.step(act)
        if info["meta"]["steps_counter"] % 3 == 0:
            # each environments only has 3 evals and will thus stop
            # even if wait period is 5
            # Will end at -1 since it does another step
            assert info["evals_left"] == 0
            assert "final_observation" in info
            assert "final_r" in info
        if info["meta"]["steps_counter"] % 3 == 1:
            assert info["evals_left"] == 2
            assert "final_observation" not in info
            assert "final_r" not in info
    assert info["meta"]["steps_counter"] == 9
    assert info["meta"]["timer"] == 9


def test_force_stop_done(ec_env):
    if not isinstance(ec_env.env_list[0].func, Sphere):
        pytest.skip("only sphere due to known optimum")

    _, info = ec_env.reset(options={"online": True})
    while not ec_env.stop():
        act = [0, 0, 0]  # This is optimal and will stop the env
        _, _, _, _, info = ec_env.step(act)
        # This costs 2 evaluations since it will auto-reset and step again
        assert info["evals_left"] == 1
        assert "final_observation" in info
        # no final_r since done on the first step
        assert "final_r" not in info
    assert info["meta"]["steps_counter"] == 3
    assert info["meta"]["timer"] == 3

    # check that even resetting does not help if the envs are done
    with pytest.raises(ValueError):
        obs, info = ec_env.reset()  # Full reset, but EC environments will throw errors
    assert ec_env.unwrapped._timer == 0
    assert ec_env.unwrapped._timer == 0
    assert ec_env.stop()
