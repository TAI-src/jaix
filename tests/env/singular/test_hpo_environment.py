from jaix.env.singular import HPOEnvironmentConfig, HPOEnvironment
from jaix.env.utils.hpo import TaskType
import pytest
from ttex.config import ConfigurableObjectFactory as COF
import json


@pytest.fixture(scope="session", autouse=True)
def skip_remaining_tests():
    try:
        import tabrepo  # noqa: F401

        assert HPOEnvironment is not None
    except ImportError:
        assert HPOEnvironment is None
        pytest.skip(
            "Skipping HPO tests. If this is unexpected, check that the tabrepo extra is installed."
        )


@pytest.fixture
def env():
    config = HPOEnvironmentConfig(
        training_budget=50,
        task_type=TaskType.C1,
        repo_name="D244_F3_C1530_30",
        cache=True,
    )
    env = COF.create(HPOEnvironment, config, func=0, inst=0)
    return env


@pytest.mark.skip("Just for checking")
def test_one_max():
    config = HPOEnvironmentConfig(
        training_budget=500,
        task_type=TaskType.C1,
        repo_name="D244_F3_C1530_30",
        cache=True,
    )
    worse = 0
    for f in HPOEnvironment.info(config)["funcs"]:
        for i in range(3):
            env = COF.create(HPOEnvironment, config, func=f, inst=i)
            env.reset(options={"online": True})

            act = env.action_space.sample()
            one_max = [1] * len(act)

            obs_act, _, _, _, _ = env.step(act)
            obs_one_max, _, _, _, _ = env.step(one_max)

            if obs_one_max[0] >= obs_act[0]:
                worse += 1
                return  # To shorten the test
    assert worse >= 0


def test_init(env):
    assert env.training_time == 0
    assert env.action_space.n == len(env.tabrepo_adapter.configs)
    assert env.training_budget == 50


def test_step(env):
    env.reset(options={"online": True})
    assert env.num_resets == 1

    obs, r, term, trunc, info = env.step([0] * env.action_space.n)
    assert obs in env.observation_space
    assert r is None
    assert not term
    assert not trunc
    assert info["env_step"] == 0

    obs, r, term, trunc, info = env.step(env.action_space.sample())
    assert obs in env.observation_space
    assert r is None
    assert not trunc
    assert info["env_step"] > 0


def test_stop(env):
    env.reset(options={"online": True})
    assert not env.stop()
    while not env.stop():
        obs, r, _, _, info = env.step(env.action_space.sample())
        ensembles = json.loads(info["ensembles"])
        assert len(ensembles[str(obs[0])]) >= 1
    assert env.training_budget <= env.training_time
    assert r is None


def test_instance_seeding():
    config = HPOEnvironmentConfig(
        training_budget=500,
        task_type=TaskType.C1,
        repo_name="D244_F3_C1530_30",
        cache=True,
    )
    env1 = COF.create(HPOEnvironment, config, func=0, inst=0)
    env2 = COF.create(HPOEnvironment, config, func=0, inst=0)
    env3 = COF.create(HPOEnvironment, config, func=0, inst=1)

    act = env1.action_space.sample()
    obs1, _, _, _, _ = env1.step(act)
    obs2, _, _, _, _ = env2.step(act)
    obs3, _, _, _, _ = env3.step(act)

    assert obs1 == obs2
    assert obs1 != obs3
