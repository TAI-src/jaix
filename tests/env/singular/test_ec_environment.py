from jaix.env.singular.ec_env import (
    ECEnvironment,
    ECEnvironmentConfig,
)
from ttex.config import ConfigurableObjectFactory as COF
from jaix.env.utils.problem.sphere import Sphere, SphereConfig
import pytest
from . import test_handler
import numpy as np


@pytest.fixture(scope="function")
def env():
    func_config = SphereConfig(
        dimension=3,
        num_objectives=2,
        mult=1,
        x_shifts=[[0, 0, 0], [0, 0, 0]],
        y_shifts=[0, 0],
        precision=1e-8,
    )
    func = COF.create(Sphere, func_config, 1)

    config = ECEnvironmentConfig(budget_multiplier=1)
    env = COF.create(ECEnvironment, config, func, 0, 1)
    return env


def test_init(env):
    assert env.num_resets == 0
    assert type(env.func) == Sphere

    info = env._get_info()
    assert info["num_resets"] == 0
    assert info["evals_left"] == 3
    assert "ECEnvironment/Sphere/0/1" in str(env)


def test_reset(env):
    obs, info = env.reset(options={"online": True})
    assert info["num_resets"] == 1
    _, _, _, _, info = env.step([0, 0, 0])
    assert info["num_resets"] == 1
    assert info["evals_left"] == 3 - 1
    obs, info = env.reset(options={"online": True})
    assert info["num_resets"] == 2
    assert info["evals_left"] == 3 - 1


def test_step_final_target_hit(env):
    obs, reward, terminated, truncated, info = env.step([0, 0, 0])
    assert obs == env.func.min_values
    assert terminated
    assert not truncated
    assert info["evals_left"] == 3 - 1
    assert env.stop() == True


def test_step_no_evals_left(env):
    obs, reward, terminated, truncated, info = env.step([1, 1, 1])
    assert not terminated
    assert not truncated
    assert info["evals_left"] == 3 - 1
    assert env.stop() == False
    obs, reward, terminated, truncated, info = env.step([1, 1, 1])
    assert not terminated
    assert not truncated
    assert info["evals_left"] == 3 - 2
    assert env.stop() == False
    obs, reward, terminated, truncated, info = env.step([1, 1, 1])
    assert not terminated
    assert truncated
    assert info["evals_left"] == 3 - 3
    assert env.stop() == True


"""
def test_step_out_of_bounds(env):
    # Deactivated truncation for EC
    obs, reward, terminated, truncated, info = env.step([-10, 0, 0])
    assert not terminated
    assert truncated
    assert info["evals_left"] == 3 - 1
    assert env.stop() == False
"""


def test_render(env):
    test_handler.last_record = None
    env.render()
    assert test_handler.last_record is not None
