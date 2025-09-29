import pytest
from ttex.config import ConfigurableObjectFactory as COF
from jaix.env.utils.problem.sphere import Sphere, SphereConfig
from jaix.suite.ec_suite import ECSuite, ECSuiteConfig
from jaix.suite.suite import AggType
from jaix.env.singular.ec_env import ECEnvironmentConfig
import os


@pytest.fixture(scope="function")
def func_config():
    func_config = SphereConfig(
        dimension=3,
        num_objectives=2,
        mult=1,
        x_shifts=[[0, 0, 0], [0, 0, 0]],
        y_shifts=[0, 0],
        precision=1e-8,
    )
    return func_config


@pytest.fixture(scope="function")
def env_config():
    config = ECEnvironmentConfig(budget_multiplier=1)
    return config


def test_init(func_config, env_config):
    print("fasf")
    config = ECSuiteConfig(
        [Sphere],
        [func_config],
        env_config,
        instances=list(range(2)),
        agg_instances=1,
    )
    suite = COF.create(ECSuite, config)

    assert suite.func_classes[0] == Sphere
    assert suite.func_configs[0].dimension == 3


def test_get_envs(func_config, env_config):
    config = ECSuiteConfig(
        [Sphere],
        [func_config],
        env_config,
        instances=list(range(1)),
        agg_instances=1,
    )
    suite = COF.create(ECSuite, config)

    inst = 0
    for env in suite.get_envs():
        assert isinstance(env.unwrapped.func, Sphere)
        assert not env.stop()
        rec_file = env.close()
        assert f"ECEnvironment/Sphere/0/{inst}" in str(env)
        if rec_file is not None:
            os.remove(rec_file)
        inst += 1


def test_get_envs_agg(func_config, env_config):
    config = ECSuiteConfig(
        [Sphere],
        [func_config],
        env_config,
        instances=list(range(3)),
        agg_instances=1,
    )
    suite = COF.create(ECSuite, config)

    for envs in suite.get_agg_envs(AggType.INST):
        assert len(envs) == 3
        assert all([isinstance(env.func, Sphere) for env in envs])
        assert all([not env.stop() for env in envs])
        rec_files = [env.close() for env in envs]
        for rec_file in rec_files:
            if rec_file is not None:
                os.remove(rec_file)
