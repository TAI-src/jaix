from jaix import EnvironmentConfig, CompositeEnvironmentConfig
from jaix.suite import AggType, ECSuite, ECSuiteConfig
from jaix.env.composite import SwitchingEnvironment, SwitchingEnvironmentConfig
from jaix.env.utils.switching_pattern import (
    SeqRegSwitchingPatternConfig,
    SeqRegSwitchingPattern,
)
from jaix.env.singular import ECEnvironmentConfig, ECEnvironment
from jaix.env.utils.problem import Sphere, SphereConfig
import pytest
from jaix import EnvironmentFactory as EF
import gymnasium as gym
from . import DummyWrapper, DummyWrapperConfig


@pytest.fixture(scope="function")
def comp_config():
    sp_config = SeqRegSwitchingPatternConfig(wait_period=3)
    config = SwitchingEnvironmentConfig(
        SeqRegSwitchingPattern, sp_config, real_time=False
    )
    comp_config = CompositeEnvironmentConfig(
        agg_type=AggType.INST,
        comp_env_class=SwitchingEnvironment,
        comp_env_config=config,
    )
    return comp_config


@pytest.fixture(scope="function")
def ec_config():
    func_config = SphereConfig(
        dimension=3,
        num_objectives=2,
        mult=1,
        x_shifts=[[0, 0, 0], [0, 0, 0]],
        y_shifts=[0, 0],
        precision=1e-8,
    )
    ec_config = ECEnvironmentConfig(budget_multiplier=1)
    config = ECSuiteConfig(
        [Sphere], [func_config], ec_config, instances=[0], agg_instances=1
    )
    return config


def env_config(ec_config, wrappers=None, comp_config=None, seed=None):
    env_config = EnvironmentConfig(
        suite_class=ECSuite,
        suite_config=ec_config,
        env_wrappers=wrappers,
        comp_config=comp_config,
        seed=seed,
    )
    return env_config


def test_wrapper_setup_teardown(ec_config):
    # wrappers = [(DummyWrapper, )]
    config = env_config(
        ec_config,
        wrappers=[(DummyWrapper, DummyWrapperConfig())],
        comp_config=None,
        seed=42,
    )
    assert config.setup()
    for _, wrapper_conf in config.env_wrappers:
        if isinstance(wrapper_conf, DummyWrapperConfig):
            assert wrapper_conf._stp
            assert not wrapper_conf.trdwn
    assert config.teardown()
    for _, wrapper_conf in config.env_wrappers:
        if isinstance(wrapper_conf, DummyWrapperConfig):
            assert wrapper_conf.trdwn


def test_singular(ec_config):
    config = env_config(ec_config)
    for env in EF.get_envs(config):
        assert isinstance(env, gym.Wrapper)
        assert isinstance(env.unwrapped, ECEnvironment)
        action = env.action_space.sample()
        env.step(action)
        env.close()


def test_composite(ec_config, comp_config):
    config = env_config(ec_config, comp_config=comp_config)
    for env in EF.get_envs(config):
        assert isinstance(env, gym.Wrapper)
        assert isinstance(env.unwrapped, SwitchingEnvironment)
        action = env.action_space.sample()
        env.step(action)
        env.close()
