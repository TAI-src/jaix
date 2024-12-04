from jaix.env.wrapper import PassthroughWrapper
import gymnasium as gym
from ttex.config import ConfigurableObject, Config
import logging


class LoggingWrapperConfig(Config):
    def __init__(
        self,
        logger_name: str,
        passthrough: bool = True,
    ):
        self.logger_name = logger_name
        self.passthrough = passthrough


class LoggingWrapper(PassthroughWrapper, ConfigurableObject):
    def __init__(self, config: LoggingWrapperConfig, env: gym.Env):
        ConfigurableObject.__init__(self, config)
        PassthroughWrapper.__init__(self, env, self.passthrough)
        self.logger = logging.getLogger(self.logger_name)
        self.steps = 0

    def reset(self, **kwargs):
        obs, info = self.env.reset(**kwargs)
        self.steps = 0
        return obs, info

    def step(self, action):
        (
            obs,
            r,
            term,
            trunc,
            info,
        ) = self.env.step(action)
        self.steps += 1
        self.logger.info({"env": str(self.env), "step": self.steps, "reward": r.item()})
        return obs, r, term, trunc, info
