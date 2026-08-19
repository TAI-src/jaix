import gymnasium as gym
from ttex.config import Config, ConfigurableObject
from ttex.config import ConfigurableObjectFactory as COF

from jaix.env.utils.archive.action_space import ArchiveActionSpace
from jaix.env.wrapper.archive_wrapper import ArchiveWrapper, ArchiveWrapperConfig


class ArchiveActionWrapperConfig(Config):
    """
    Configuration for the ArchiveActionWrapper.
    """

    def __init__(
        self,
        archive_wrapper_config: ArchiveWrapperConfig,
        action_space_class: type[ArchiveActionSpace],
        action_space_config: Config,
    ):
        super().__init__()
        self.archive_wrapper_config = archive_wrapper_config
        self.action_space_class = action_space_class
        self.action_space_config = action_space_config


class ArchiveActionWrapper(ArchiveWrapper):
    """
    A wrapper that modifies the action space to include actions on a population.
    """

    config_class = ArchiveActionWrapperConfig  # type: ignore[assignment]

    def __init__(self, config: ArchiveActionWrapperConfig, env: gym.Env, **kwargs):
        ConfigurableObject.__init__(self, config)
        ArchiveWrapper.__init__(
            self, self.archive_wrapper_config, env, **kwargs
        )  # This creates the archive and sets self.archive
        self.action_space = COF.create(
            self.action_space_class, self.action_space_config, archive=self.archive
        )

    def step(self, action):
        # Translate the action into the archive's action space
        translated_action = self.action_space.translate(action)
        obs, reward, term, trunc, info = super().step(translated_action)
        info["archive_action"] = translated_action
        info["original_action"] = action
        return obs, reward, term, trunc, info
