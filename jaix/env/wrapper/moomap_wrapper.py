from jaix.env.wrapper.archive_wrapper import ArchiveWrapper, ArchiveWrapperConfig
from jaix.env.utils.mo_sizing import get_num_refpoints, get_ref_dirs
from ttex.config import Config, ConfigurableObject


class MoomapWrapperConfig(Config):
    """
    Configuration for the MoomapWrapper.
    """

    def __init__(
        self,
        archive_class: type[ArchiveWrapper],
        archive_config: Config | None,
        replace_reward: bool = True,
        passthrough: bool = True,
        num_refpoints: int | str = "original",
    ):
        super().__init__(
            archive_class=archive_class,
            archive_config=archive_config,
            replace_reward=replace_reward,
            passthrough=passthrough,
        )
        if isinstance(num_refpoints, int):
            assert num_refpoints > 0, "num_refpoints must be a positive integer"
            self.num_refpoints = num_refpoints
        elif isinstance(num_refpoints, str):
            self.num_refpoints = -1


class MoomapWrapper(ArchiveWrapper, ConfigurableObject):
    """
    A wrapper that modifies the action space to include actions on a population.
    """

    config_class = MoomapWrapperConfig

    def __init__(self, config: ArchiveWrapperConfig, env, **kwargs):
        ArchiveWrapper.__init__(self, config, env, **kwargs)
