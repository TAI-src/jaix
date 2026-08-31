from abc import ABC, abstractmethod
from typing import Any

import numpy as np
from jaix.env.wrapper.archive_wrapper import EnvironmentStepEntry
from moocore import pareto_rank, is_nondominated, hv_approx, hypervolume
from ttex.config import Config, ConfigurableObject

from jaix.env.utils.archive.archive import Archive, ArchiveEntry
import gymnasium as gym
from jaix.env.singular.ec_env import get_ideal_nadir
from enum import Enum


class MOArchiveEntry(EnvironmentStepEntry[np.ndarray], ABC):
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
        self.hv_contribution: float = np.nan
        self.rank: int = -1

    def parse(self) -> np.ndarray:
        return self.obs

    @property
    def objectives(self) -> np.ndarray:
        return self.parse()


class SecondaryCriterion(ABC):
    @abstractmethod
    def rank(self, entries: list[MOArchiveEntry]) -> list[int]: ...


class KeepDominated(Enum):
    ALL = "all"
    NONE = "none"
    STRICT = "strict"


class MOArchiveConfig(Config):
    def __init__(
        self,
        archive_entry_class: type[MOArchiveEntry],
        secondary_criterion: SecondaryCriterion,
        max_size: int | None = None,
        keep_dominated: KeepDominated = KeepDominated.NONE,
        only_new_entries: bool = False,
        hv_approx_samples: int | None = 262_144,
    ) -> None:
        super().__init__()
        self.max_size = max_size
        self.archive_entry_class = archive_entry_class
        self.secondary_criterion = secondary_criterion
        self.keep_dominated = keep_dominated
        self.only_new_entries = only_new_entries
        self.hv_approx_samples = hv_approx_samples


class MOArchive(Archive, ConfigurableObject):
    config_class = MOArchiveConfig

    def __init__(self, config: MOArchiveConfig, env: gym.Env, **kwargs) -> None:
        ConfigurableObject.__init__(self, config)
        Archive.__init__(self, max_size=config.max_size)
        self.archived_entries: list[MOArchiveEntry] = []

        self.ideal_point, self.nadir_point = get_ideal_nadir(env)

    def _add(self, entries: list[ArchiveEntry]) -> list[dict[str, Any]]:
        """
        Add an entry to the archive. This method is meant to be implemented by subclasses.
        """
        assert all(
            isinstance(entry, MOArchiveEntry) for entry in entries
        ), "All entries must be of type MOArchiveEntry"
        if self.only_new_entries:
            # Only new entries are added to the archive, so we replace the archived entries with the new ones
            self.archived_entries = entries
        else:
            # All entries are added to the archive, so we extend the archived entries with the new ones
            self.archived_entries.extend(entries)

        if self.keep_dominated != KeepDominated.ALL:
            # Only keep non-dominated entries in the archive
            MOArchive.remove_dominated(self.archived_entries, self.keep_dominated)

        if self.max_size is not None and len(self.archived_entries) > self.max_size:
            # We need to remove some entries from the archive to keep it within the max size
            max_rank = MOArchive.crit_rank(self.archived_entries, self.max_size)
            self.archived_entries = [
                entry for entry in self.archived_entries if entry.rank < max_rank
            ]
            critical_entries = [
                entry for entry in self.archived_entries if entry.rank == max_rank
            ]
            # TODO: SecondaryCriterion

        raise NotImplementedError("Subclasses must implement the _add method.")

    @property
    def archive_entry_type(self) -> type[MOArchiveEntry]:
        return self.config.archive_entry_class

    @property
    def size(self) -> int:
        return len(self.archived_entries)

    @staticmethod
    def remove_dominated(
        entries: list[MOArchiveEntry], keep_dominated: KeepDominated
    ) -> list[MOArchiveEntry]:
        """
        Remove dominated entries from the archive.
        """
        points = np.array([entry.objectives for entry in entries])
        keep_weakly = keep_dominated != KeepDominated.STRICT
        non_dominated_mask = is_nondominated(points, keep_weakly=keep_weakly)
        entries = [entry for entry, keep in zip(entries, non_dominated_mask) if keep]
        return entries

    @staticmethod
    def crit_rank(entries: list[MOArchiveEntry], max_size: int) -> int:
        """
        Returns the critical rank of the archive given the current entries and the maximum size.
        Also adds the rank to each entry in the archive.
        """
        points = np.array([entry.objectives for entry in entries])
        ranks = pareto_rank(points)
        for entry in entries:
            entry.rank = ranks[entries.index(entry)]
        rank = 0
        added = 0
        while added < max_size and rank <= max(ranks):
            size_rank = sum(ranks == rank)
            if added + size_rank <= max_size:
                rank += 1
                added += size_rank
            else:
                break
        return rank

    @property
    def score(self) -> float:
        """
        Return the score of the archive as a float.
        The score is the hypervolume of the non-dominated solutions in the archive.
        """
        if self.hv_approx_samples is None:
            hv = hypervolume(self.archived_entries, ref=self.nadir_point)
        else:
            hv = hv_approx(self.archived_entries, ref=self.nadir_point)
        return hv

    def get_archive_stats(self) -> dict[str, Any]:
        """
        Return a dictionary with the current archive stats.
        """
        return {
            "size": self.size,
            "score": self.score,
        }
