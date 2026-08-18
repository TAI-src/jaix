import gymnasium as gym
from ttex.config import Config, ConfigurableObject
from ttex.config import ConfigurableObjectFactory as COF

from jaix.env.utils.archive.archive import Archive
from jaix.env.wrapper.passthrough_wrapper import PassthroughWrapper


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
    A wrapper that maintains an archive of solutions based on the observations and rewards from the environment.
    This is just for recording purposes, and does not affect the environment's behavior.
    The reward is based on the archive's evaluation of the solution, which may be different from the environment's reward.
    """

    config_class = ArchiveWrapperConfig

    def __init__(self, config: ArchiveWrapperConfig, env: gym.Env, **kwargs):
        ConfigurableObject.__init__(self, config)
        PassthroughWrapper.__init__(self, env, config.passthrough)
        if config.archive_config is None:
            self.archive = config.archive_class(**kwargs)
        else:
            self.archive = COF.create(
                config.archive_class, config.archive_config, **kwargs
            )

    def reset(self, **kwargs):
        self.archive.reset()
        # TODO: pop the used kwargs and pass the rest on
        return self.env.reset(**kwargs)

    def step(self, action):
        obs, r, term, trunc, info = self.env.step(action)
        # We assume that the action is a solution to be added to the archive, and the reward is the fitness of that solution.
        archive_reward = self.archive.add(action, r)
        ret_reward = archive_reward if self.replace_reward else r
        return obs, ret_reward, term, trunc, info
