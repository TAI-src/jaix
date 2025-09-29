from jaix.suite.suite import Suite, SuiteConfig, AggType
from . import DummyEnvConfig, DummyConfEnv
from collections.abc import Iterable
import pytest


def init_suite():
    config = SuiteConfig(
        env_class=DummyConfEnv,
        env_config=DummyEnvConfig(dimension=6),
        functions=[0, 1],
        instances=list(range(5)),
        agg_instances=[0, 1, 2],
    )
    return Suite(config)


@pytest.mark.parametrize("inst,comp_env_num", [(20, 2), (3, 3)])
def test_suite_config_agg_instances_many(inst, comp_env_num):
    config_dict = {
        "env_class": DummyConfEnv,
        "env_config": DummyEnvConfig(dimension=6),
        "functions": [0, 1],
        "instances": list(range(inst)),
        "comp_env_num": comp_env_num,
    }
    if inst == 20:
        agg_options = [
            3,  # 3 random permutations
            [0, 17, 25],  # specific indices
            [(0, 1), (2, 3), (4, 5)],  # specific permutations
        ]
    else:
        agg_options = [
            6,  # all random permutations
            list(range(6)),  # specific indices
            None,  # all permutations
        ]
    for agg in agg_options:
        config = SuiteConfig(**config_dict, agg_instances=agg)
        assert isinstance(
            config.agg_instances, list
        ), "agg_instances should be iterable"
        assert len(config.agg_instances) == agg_options[0]
        assert all(
            len(i) == comp_env_num for i in config.agg_instances
        ), "agg_instances should contain comp_env_num instances"
        assert all(
            all(isinstance(x, int) for x in i) for i in config.agg_instances
        ), "agg_instances should contain only integers"
        assert len(set(config.agg_instances)) == len(
            config.agg_instances
        ), "agg_instances should contain unique tuples"


def test_init():
    suite = init_suite()
    counter = 0
    assert len(suite.agg_instances) == 3


def test_get_envs():
    suite = init_suite()
    func = 0
    inst = 0
    for env in suite.get_envs():
        assert isinstance(env, DummyConfEnv)
        assert not env.stop()
        env.step(env.action_space.sample())
        env.close()
        assert env.inst == inst
        assert env.func == func
        inst += 1
        if inst == len(suite.instances):
            inst = 0
            func += 1


def test_get_agg_envs():
    suite = init_suite()
    counter = 0
    for envs in suite.get_agg_envs(AggType.INST, seed=5):
        assert len(envs) == len(suite.instances)
        assert isinstance(envs[0], DummyConfEnv)
        counter += 1
    assert counter == len(suite.functions) * len(
        suite.agg_instances
    ), "Number of aggregated environments does not match expected count"
