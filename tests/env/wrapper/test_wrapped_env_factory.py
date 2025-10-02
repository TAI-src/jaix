from jaix.env.wrapper.wrapped_env_factory import WrappedEnvFactory as WEF
from gymnasium.wrappers import RescaleAction
from . import DummyEnv
import gymnasium as gym
from ttex.config import ConfigurableObject, Config


def test_plain_wrap():
    base_env = DummyEnv()
    assert base_env.action_space.low[0] == -5
    assert base_env.action_space.high[0] == 5

    wrappers = [(RescaleAction, {"min_action": 0, "max_action": 1})]

    wrapped_env = WEF.wrap(base_env, wrappers)
    assert wrapped_env.action_space.low[0] == 0
    assert wrapped_env.action_space.high[0] == 1


class DummyWrapperConfig(Config):
    def __init__(self, test_value):
        self.test_value = test_value


class DummyWrapper(gym.Wrapper, ConfigurableObject):
    config_class = DummyWrapperConfig

    def __init__(self, config: DummyWrapperConfig, env: gym.Env):
        ConfigurableObject.__init__(self, config)
        gym.Wrapper.__init__(self, env)


def test_config_wrap():
    base_env = DummyEnv()
    assert not hasattr(base_env, "test_value")

    wrappers = [(DummyWrapper, DummyWrapperConfig(test_value=7))]
    wrapped_env = WEF.wrap(base_env, wrappers)
    assert wrapped_env.test_value == 7


def test_config_wrap_dict():
    base_env = DummyEnv()
    assert not hasattr(base_env, "test_value")

    wrappers = [
        (
            DummyWrapper,
            {
                "tests.env.wrapper.test_wrapped_env_factory.DummyWrapperConfig": {
                    "test_value": 9
                }
            },
        )
    ]
    wrapped_env = WEF.wrap(base_env, wrappers)
    # Check that the config was converted to object
    assert isinstance(wrappers[0][1], DummyWrapperConfig)
    # Check that the env was wrapped correctly
    assert wrapped_env.test_value == 9
