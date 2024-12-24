from jaix.env.singular import HPOEnvironmentConfig, HPOEnvironment
from jaix.env.utils.hpo import TaskType
from ttex.config import ConfigurableObjectFactory as COF
import pytest


@pytest.fixture
def env():
    config = HPOEnvironmentConfig(
        training_budget=10000,
        task_type=TaskType.C1,
        repo_name="D244_F3_C1530_3",
        load_predictions=False,
        cache=True,
    )
    env = COF.create(HPOEnvironment, config, 0)
    return env


def test_init(env):
    assert env.training_time == 0
    assert env.training_budget == 10000
    assert len(env.action_space.nvec) == 1
    assert env.action_space.nvec[0] > 0


def test_step(env):
    env.reset(options={"online": True})
    assert env.num_resets == 1

    obs, r, term, trunc, info = env.step(env.action_space.sample())
    assert obs in env.observation_space
    assert r == obs[0]
    assert not term
    assert not trunc
    assert info["training_time"] > 0


def test_stop(env):
    env.reset(options={"online": True})
    assert not env.stop()
    while not env.stop():
        env.step(env.action_space.sample())
    assert env.training_budget <= env.training_time
