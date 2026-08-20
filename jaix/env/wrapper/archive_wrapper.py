from abc import ABC
from typing import Any, TypeVar

import gymnasium as gym
from ttex.config import Config, ConfigurableObject
from ttex.config import ConfigurableObjectFactory as COF

from jaix.env.utils.archive.archive import Archive, ArchiveEntry
from jaix.env.wrapper.passthrough_wrapper import PassthroughWrapper

T = TypeVar("T")


class EnvironmentStepEntry(ArchiveEntry[T], ABC):
    def __init__(
        self,
        action: Any,
        obs: Any,
        reward: float,
        terminated: bool,
        truncated: bool,
        info: dict,
    ):
        self.action = action
        self.obs = obs
        self.reward = reward
        self.terminated = terminated
        self.truncated = truncated
        self.info = info


class ArchiveWrapperConfig(Config):
    def __init__(
        self,
        archive_class: type[Archive],
        archive_config: Config | None,
        replace_reward: bool = True,
        passthrough: bool = True,
    ):
        Config.__init__(self)
        self.archive_class = archive_class
        self.archive_config = archive_config
        self.replace_reward = replace_reward
        self.passthrough = passthrough


class ArchiveWrapper(ConfigurableObject, PassthroughWrapper):
    """
    A wrapper that maintains an archive with entries that record the environment's step information.
    This is just for recording purposes, and does not affect the environment's behavior.
    The reward is based on the archive's evaluation of the solution, which may be different from the environment's reward.
    """

    config_class = ArchiveWrapperConfig

    def __init__(self, config: ArchiveWrapperConfig, env: gym.Env, **kwargs):
        ConfigurableObject.__init__(self, config)
        PassthroughWrapper.__init__(self, env, config.passthrough)
        kwargs.setdefault("env", env)
        if config.archive_config is None:
            self.archive = config.archive_class(**kwargs)
        else:
            self.archive = COF.create(
                config.archive_class, config.archive_config, **kwargs
            )
        assert issubclass(self.archive.archive_entry_type, EnvironmentStepEntry)
        self.environment_step_class: type[EnvironmentStepEntry] = (
            self.archive.archive_entry_type
        )

    def reset(self, **kwargs):
        self.archive.reset()
        return self.env.reset(**kwargs)

    def step(self, action):
        obs, r, term, trunc, info = self.env.step(action)
        # We assume that the action is a solution to be added to the archive, and the reward is the fitness of that solution.
        entry = self.environment_step_class(
            action=action,
            obs=obs,
            reward=r,
            terminated=term,
            truncated=trunc,
            info=info,
        )
        archive_reward = self.archive.add(entry)
        ret_reward = archive_reward if self.replace_reward else r
        return obs, ret_reward, term, trunc, info
