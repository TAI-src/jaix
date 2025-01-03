from jaix.suite import Suite, SuiteConfig, AggType
from . import DummyEnvConfig, DummyConfEnv


def init_suite():
    config = SuiteConfig(
        env_class=DummyConfEnv,
        env_config=DummyEnvConfig(dimension=6),
        num_instances=5,
        num_agg_instances=3,
    )
    return Suite(config)


def test_init():
    suite = init_suite()
    assert suite.num_instances == 5


def test_get_envs():
    suite = init_suite()
    counter = 0
    for env in suite.get_envs():
        assert isinstance(env, DummyConfEnv)
        assert not env.stop()
        env.step(env.action_space.sample())
        env.close()
        assert env.inst == counter
        counter += 1
    assert counter == suite.num_instances


def test_get_agg_envs():
    suite = init_suite()
    counter = 0
    insts = []
    for envs in suite.get_agg_envs(AggType.INST, seed=5):
        assert len(envs) == suite.num_agg_instances
        assert isinstance(envs[0], DummyConfEnv)
        insts.append(envs[0].inst)
        counter += 1
    assert counter == suite.num_instances
    assert not all([inst == insts[0] for inst in insts])
