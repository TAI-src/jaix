from jaix.env.wrapper.passthrough_wrapper import PassthroughWrapper
from ttex.config import Config, ConfigurableObject, ConfigurableObjectFactory as COF
from jaix.env.archive.archive import Archive
import gymnasium as gym


class ArchiveActionWrapperConfig(Config):
    """
    Configuration for the PopulationActionWrapper.
    """

    def __init__(
        self,
        archive_class: type[Archive],
        archive_config: Config,
        passthrough: bool = True,
    ):
        super().__init__()
        self.archive_class = archive_class
        self.archive_config = archive_config
        self.passthrough = passthrough


class ArchiveActionWrapper(ConfigurableObject, PassthroughWrapper):
    """
    A wrapper that modifies the action space to include actions on a population.
    """

    config_class = ArchiveActionWrapperConfig

    def __init__(self, config: ArchiveActionWrapperConfig, env: gym.Env, **kwargs):
        ConfigurableObject.__init__(self, config)
        PassthroughWrapper.__init__(self, env, passthrough=self.passthrough)
        self.archive = COF.create(self.archive_class, self.archive_config, **kwargs)
