from jaix.env.wrapper.passthrough_wrapper import PassthroughWrapper
import gymnasium as gym
from ttex.config import Config, ConfigurableObject


class DummyWrapperConfig(Config):
    def __init__(self, passthrough: bool = True):
        self.passthrough = passthrough
        self._stp = False
        self.trdwn = False

    def _setup(self):
        self._stp = True
        return self._stp

    def _teardown(self):
        self.trdwn = True
        return self.trdwn


class DummyWrapper(PassthroughWrapper, ConfigurableObject):
    config_class = DummyWrapperConfig

    def __init__(self, config: Config, env: gym.Env):
        ConfigurableObject.__init__(self, config)
        PassthroughWrapper.__init__(self, env, passthrough=self.passthrough)
        self.stop_dict = {}
        self.env_steps = 0

    def _stop(self):
        return self.stop_dict

    def step(self, action):
        obs, r, term, trunc, info = self.env.step(action)
        self.env_steps += 1
        info["env_step"] = self.env_steps + 1
        return obs, r, term, trunc, info
