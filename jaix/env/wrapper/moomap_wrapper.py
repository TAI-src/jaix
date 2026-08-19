from jaix.env.wrapper.archive_wrapper import ArchiveWrapper, ArchiveWrapperConfig
from jaix.env.utils.mo_sizing import get_num_refpoints, get_ref_dirs
from ttex.config import Config, ConfigurableObject
from jaix.env.utils.archive.bin_archive import BinArchive, BinArchiveConfig
from jaix.env.utils.archive.rv_binning_strategy import (
    RVBinningStrategy,
    RVBinningStrategyConfig,
)
from jaix.env.singular.ec_env import ECEnvironment
import gymnasium as gym
from jaix.env.utils.problem.static_problem import StaticFunction


class MoomapArchiveConfig(Config):
    def __init__(
        self,
        np_bin: int,
        coverage_weight: float,
        allow_close_elites: bool = True,
        num_refpoints: int | str = "original",
    ):
        super().__init__()
        self.np_bin = np_bin
        self.coverage_weight = coverage_weight
        self.allow_close_elites = allow_close_elites
        self.num_refpoints = num_refpoints


class MoomapArchiveEntry:
    """
    An entry in the MoomapArchive.
    """

    # TODO: observation and distance to ideal as the fitness


class MoomapArchive(BinArchive):
    """
    A BinArchive that uses reference directions for binning.
    """

    def __init__(self, config: MoomapArchiveConfig, func: StaticFunction):
        assert func.num_objectives > 1, "Function must be multi-objective"
        assert hasattr(func, "nadir_point"), "Function must have nadir_point attribute"
        assert hasattr(func, "ideal_point"), "Function must have ideal_point attribute"

        # Determine the number of reference points
        if isinstance(config.num_refpoints, int):
            num_refpoints = config.num_refpoints
        elif isinstance(config.num_refpoints, str):
            num_refpoints = get_num_refpoints(func.num_objectives)
        else:
            raise ValueError(
                "num_refpoints must be a positive integer or a known string"
            )

        # Generate reference directions
        ref_dirs = get_ref_dirs(func.num_objectives, num_refpoints)

        # Create the binning strategy
        binning_strategy_config = RVBinningStrategyConfig()
        binning_strategy = RVBinningStrategy(
            binning_strategy_config, ref_dirs, func.ideal_point, func.nadir_point
        )

        # Initialize the BinArchive with the binning strategy
        super().__init__(
            binning_strategy,
            config.np_bin,
            config.coverage_weight,
            config.allow_close_elites,
        )


class MoomapWrapperConfig(Config):
    """
    Configuration for the MoomapWrapper.
    """

    def __init__(
        self,
        np_bin: int,
        coverage_weight: float,
        allow_close_elites: bool = True,
        replace_reward: bool = True,
        passthrough: bool = True,
        num_refpoints: int | str = "original",
    ):
        super().__init__()
        self.archive_class = BinArchive
        self.binning_strategy_config = RVBinningStrategyConfig
        self.binning_strategy_class = RVBinningStrategy

        if isinstance(num_refpoints, int):
            assert num_refpoints > 0, "num_refpoints must be a positive integer"
            self.num_refpoints = num_refpoints
        elif isinstance(num_refpoints, str):
            self.num_refpoints = -1


class MoomapWrapper(ConfigurableObject):
    """
    A wrapper that modifies the action space to include actions on a population.
    """

    config_class = MoomapWrapperConfig

    def __init__(self, config: ArchiveWrapperConfig, env: gym.Env):
        ConfigurableObject.__init__(self, config)
        assert isinstance(
            env.unwrapped, ECEnvironment
        ), "MoomapWrapper can only be used with ECEnvironment"
        assert hasattr(
            env.unwrapped.func, "nadir_point"
        ), "Function must have nadir_point attribute"
        assert hasattr(
            env.unwrapped.func, "ideal_point"
        ), "Function must have ideal_point attribute"
        assert env.unwrapped.func.num_objectives > 1, "Function must be multi-objective"
