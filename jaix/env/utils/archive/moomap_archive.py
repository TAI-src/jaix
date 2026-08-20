from jaix.env.singular.ec_env import ECEnvironment
from jaix.env.wrapper.archive_wrapper import EnvironmentStepEntry
from ttex.config import Config, ConfigurableObject, ConfigurableObjectFactory as COF
from jaix.env.utils.archive.bin_archive import BinArchive, BinArchiveEntry
from jaix.env.utils.archive.archive import Archive, ArchiveEntry
from jaix.env.utils.mo_sizing import get_ref_dirs
from jaix.env.utils.problem.static_problem import StaticProblem
from jaix.env.utils.archive.rv_binning_strategy import (
    RVBinningStrategy,
    RVBinningStrategyConfig,
)
import numpy as np
import gymnasium as gym


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


class MoomapArchiveEntry(
    EnvironmentStepEntry[tuple[np.ndarray, float]], BinArchiveEntry[np.ndarray]
):
    """
    An entry in the MoomapArchive.
    """

    def __init__(
        self,
        action: np.ndarray,
        obs: np.ndarray | list,
        reward: float,
        terminated: bool,
        truncated: bool,
        info: dict,
    ):
        EnvironmentStepEntry.__init__(
            self, action, np.asarray(obs), reward, terminated, truncated, info
        )
        self._fitness = (
            np.nan
        )  # This is the distance to ideal which is not known at the time of creation. It will be set when the entry is added to the archive.

    def parse(self) -> tuple[np.ndarray, float]:
        # Observation in EC environment is by convention the objective vector
        # This is used as the bin_sample for the RVBinningStrategy (reference vector based)
        return self.obs, self._fitness

    def set_fitness(self, ideal_point: np.ndarray):
        """
        Set the fitness of the entry based on the distance to the ideal point.
        """
        self._fitness = np.linalg.norm(self.obs - ideal_point)


class MoomapArchive(BinArchive):
    """
    A BinArchive that uses reference directions for binning.
    """

    config_class = MoomapArchiveConfig

    def __init__(self, config: MoomapArchiveConfig, env: gym.Env):
        assert isinstance(
            env.unwrapped, ECEnvironment
        ), "MoomapArchive can only be used with ECEnvironment"
        self.func = env.unwrapped.func
        ConfigurableObject.__init__(self, config)
        assert self.func.num_objectives > 1, "Function must be multi-objective"

        # Generate reference directions
        self.ref_dirs = get_ref_dirs(self.func.num_objectives, self.num_refpoints)
        self.num_refpoints = len(self.ref_dirs)
        self.ideal_point = getattr(self.func, "ideal_point", None)
        self.nadir_point = getattr(self.func, "nadir_point", None)
        assert self.ideal_point is not None, "Function must have ideal_point attribute"
        assert self.nadir_point is not None, "Function must have nadir_point attribute"

        # Create the binning strategy
        self.binning_config = RVBinningStrategyConfig()
        self.binning_strategy = RVBinningStrategy
        self.binner = COF.create(
            self.binning_strategy,
            config=self.binning_config,
            ref_dirs=self.ref_dirs,
            ideal=self.ideal_point,
            nadir=self.nadir_point,
        )
        # Manually set the values that are used in the BinArchive
        self.n_bins = self.num_refpoints
        self.archive_entry_class = MoomapArchiveEntry
        self.max_fitness = np.linalg.norm(self.nadir_point - self.ideal_point)
        Archive.__init__(self, max_size=self.n_bins)
        super().reset()  # Reset the BinArchive, which will also init the bins and the stats

    def add(self, entry: ArchiveEntry):
        assert isinstance(
            entry, MoomapArchiveEntry
        ), "Entry must be a MoomapArchiveEntry"
        # Set the fitness of the entry based on the distance to the ideal point
        entry.set_fitness(self.ideal_point)
        return super().add(entry)
