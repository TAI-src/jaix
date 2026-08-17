from jaix.env.wrapper.passthrough_wrapper import PassthroughWrapper
from ttex.config import Config, ConfigurableObject, ConfigurableObjectFactory as COF
from jaix.env.utils.archive.archive import Archive
import gymnasium as gym


class ArchiveWrapperConfig(Config):
    def __init__(
        self,
        archive_class: type[Archive],
        archive_config: Config,
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
        self.archive = COF.create(config.archive_class, config.archive_config, **kwargs)

    def reset(self, **kwargs):
        self.archive = COF.create(self.archive_class, self.archive_config, **kwargs)
        # TODO: pop the used kwargs and pass the rest on
        return self.env.reset()

    def step(self, action):
        obs, r, term, trunc, info = self.env.step(action)
        # Obs is assumed to be the multi-objective reward, or else some sort of descriptor of the solution. The archive will then store solutions based on the descriptor and archive
        # The reward is assumed to be the scalarized reward, or some sort of fitness value. How that will be incorporated into the archive is up to the archive implementation.
        archive_reward = self.archive.add(obs, r)
        ret_reward = archive_reward if self.replace_reward else r
        return obs, ret_reward, term, trunc, info
