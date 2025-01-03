import pytest
from ttex.config import ConfigurableObjectFactory as COF
from jaix.env.utils.problem import Sphere, SphereConfig
from jaix.suite import ECSuite, ECSuiteConfig, AggType
from jaix.env.singular import ECEnvironmentConfig
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
    config = ECSuiteConfig(
        Sphere, func_config, env_config, num_instances=1, num_agg_instances=1
    )
    suite = COF.create(ECSuite, config)

    assert suite.func_class == Sphere
    assert suite.func_config.dimension == 3


def test_get_envs(func_config, env_config):
    config = ECSuiteConfig(
        Sphere, func_config, env_config, num_instances=1, num_agg_instances=1
    )
    suite = COF.create(ECSuite, config)

    for env in suite.get_envs():
        assert isinstance(env.unwrapped.func, Sphere)
        assert not env.stop()
        rec_file = env.close()
        if rec_file is not None:
            os.remove(rec_file)


def test_get_envs_agg(func_config, env_config):
    config = ECSuiteConfig(
        Sphere, func_config, env_config, num_instances=1, num_agg_instances=3
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
