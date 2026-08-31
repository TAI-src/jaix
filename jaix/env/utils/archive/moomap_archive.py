from typing import cast

import gymnasium as gym
import numpy as np
from ttex.config import Config, ConfigurableObject
from ttex.config import ConfigurableObjectFactory as COF

from jaix.env.singular.ec_env import ECEnvironment
from jaix.env.utils.archive.archive import Archive, ArchiveEntry
from jaix.env.utils.archive.bin_archive import (
    BinArchive,
    BinArchiveEntry,
    EliteSelectionStrategy,
)
from jaix.env.utils.archive.rv_binning_strategy import (
    RVBinningStrategy,
    RVBinningStrategyConfig,
)
from jaix.env.utils.mo_sizing import get_ref_dirs
from jaix.env.wrapper.archive_wrapper import EnvironmentStepEntry


class MoomapArchiveConfig(Config):
    def __init__(
        self,
        np_bin: int | list[int],
        coverage_weight: float,
        allow_close_elites: bool = True,
        num_refpoints: int | str = "original",
        elite_selection_strategy: EliteSelectionStrategy = EliteSelectionStrategy.BEST,
    ):
        super().__init__()
        self.np_bin = np_bin
        self.coverage_weight = coverage_weight
        self.allow_close_elites = allow_close_elites
        self.num_refpoints: int | str = num_refpoints
        self.elite_selection_strategy: EliteSelectionStrategy = elite_selection_strategy


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
        self._fitness: float = (
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
        self._fitness = float(np.linalg.norm(self.obs - ideal_point))


class MoomapArchive(BinArchive):
    """
    A BinArchive that uses reference directions for binning.
    """

    config_class: type[Config] = MoomapArchiveConfig

    def __init__(self, config: MoomapArchiveConfig, env: gym.Env, **kwargs):
        assert isinstance(
            env.unwrapped, ECEnvironment
        ), "MoomapArchive can only be used with ECEnvironment"
        self.func = env.unwrapped.func
        ConfigurableObject.__init__(self, config)
        assert self.func.num_objectives > 1, "Function must be multi-objective"

        # Generate reference directions
        self.ref_dirs = get_ref_dirs(self.func.num_objectives, config.num_refpoints)
        self.num_refpoints: int = len(self.ref_dirs)
        ideal_point = getattr(self.func, "ideal_point", None)
        nadir_point = getattr(self.func, "nadir_point", None)
        assert ideal_point is not None, "Function must have ideal_point attribute"
        assert nadir_point is not None, "Function must have nadir_point attribute"
        self.ideal_point: np.ndarray = cast(np.ndarray, ideal_point)
        self.nadir_point: np.ndarray = cast(np.ndarray, nadir_point)

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
        if isinstance(self.np_bin, int):
            self.np_bin: list[int] = [self.np_bin] * self.n_bins  # type: ignore
        self.archive_entry_class = MoomapArchiveEntry
        self.max_fitness = np.linalg.norm(self.nadir_point - self.ideal_point)
        Archive.__init__(self, max_size=self.n_bins)
        super().reset()  # Reset the BinArchive, which will also init the bins and the stats

    def add(self, entries: list[ArchiveEntry]) -> float:
        for entry in entries:
            assert isinstance(
                entry, MoomapArchiveEntry
            ), "Entry must be a MoomapArchiveEntry"
            # Set the fitness of the entry based on the distance to the ideal point
            entry.set_fitness(self.ideal_point)
        return super().add(entries)
