from ..utils.hpo.test_tabrepo_adapter import repo
from jaix.env.singular import HPOEnvironmentConfig, HPOEnvironment
from jaix.env.utils.hpo import TaskType
from ttex.config import ConfigurableObjectFactory as COF
import pytest


@pytest.fixture
def env(repo):
    config = HPOEnvironmentConfig(training_budget=10000)
    env = COF.create(HPOEnvironment, config, repo, TaskType.C1, 0)
    return env


def test_init(env):
    assert env.training_time == 0
    assert env.training_budget == 10000
    assert env.action_space.n > 0


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
