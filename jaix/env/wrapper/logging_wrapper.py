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
    config_class = LoggingWrapperConfig

    def __init__(self, config: LoggingWrapperConfig, env: gym.Env):
        ConfigurableObject.__init__(self, config)
        PassthroughWrapper.__init__(self, env, self.passthrough)
        self.logger = logging.getLogger(self.logger_name)
        self.log_resets = 0
        self.log_env_steps = 0
        self.log_renv_steps = 0

    def reset(self, **kwargs):
        obs, info = self.env.reset(**kwargs)
        self.log_resets += 1
        self.log_renv_steps = 0
        return obs, info

    def step(self, action):
        (
            obs,
            r,
            term,
            trunc,
            info,
        ) = self.env.step(action)
        self.log_env_steps += 1
        self.log_renv_steps += 1
        # Log per reset
        self.logger.info(
            {
                f"env/r/{self.env}": r.item(),
                f"env/resets/{self.env}": self.log_resets,
                # f"restarts/r/{self.dim}/{self.env}/{self.log_resets}": r.item(),
                "env/step": self.log_env_steps,
                # "restarts/step": self.log_renv_steps,
            }
        )
        # TODO: Figure out what info would be helpful from all the sub-wrappers etc
        return obs, r, term, trunc, info
