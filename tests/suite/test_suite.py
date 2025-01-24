from jaix.suite import Suite, SuiteConfig, AggType
from . import DummyEnvConfig, DummyConfEnv


def init_suite():
    config = SuiteConfig(
        env_class=DummyConfEnv,
        env_config=DummyEnvConfig(dimension=6),
        functions=[0, 1],
        instances=list(range(5)),
        agg_instances=[0, 1, 2],
    )
    return Suite(config)


def test_init():
    suite = init_suite()
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
    assert counter == len(suite.functions) * len(suite.agg_instances)
