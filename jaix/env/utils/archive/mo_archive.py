from abc import ABC, abstractmethod
from enum import Enum
from typing import Any, cast

import gymnasium as gym
import numpy as np
from moocore import (
    hv_approx,
    hypervolume,
    is_nondominated,
    pareto_rank,
)
from ttex.config import Config, ConfigurableObject

from jaix.env.singular.ec_env import get_ideal_nadir
from jaix.env.utils.archive.archive import Archive, ArchiveEntry
from jaix.env.utils.archive.entry_scorer import EntryScorer
from jaix.env.utils.mo_sizing import get_ref_dirs
from jaix.utils.dirty_list import DirtyList

from pymoo.algorithms.moo.nsga3 import associate_to_niches, calc_niche_count, niching


class MOArchiveEntry(ArchiveEntry[np.ndarray], ABC):
    rank: int = -1
    secondary_score: float = np.nan
    niche: int = -1
    dist_to_ref: float = np.nan

    @abstractmethod
    def parse(self) -> np.ndarray: ...

    @property
    def objectives(self) -> np.ndarray:
        return self.parse()


class KeepDominated(Enum):
    ALL = "all"
    NONE = "none"
    STRICT = "strict"


class MOArchiveConfig(Config):
    def __init__(
        self,
        archive_entry_class: type[MOArchiveEntry],
        secondary_criterion_class: type[EntryScorer],
        max_size: int | None = None,
        keep_dominated: KeepDominated = KeepDominated.NONE,
        only_new_entries: bool = False,
        hv_approx_samples: int | None = 262_144,
        num_refpoints: int | str = "original",
    ) -> None:
        super().__init__()
        self.max_size = max_size
        self.archive_entry_class = archive_entry_class
        self.secondary_criterion_class = secondary_criterion_class
        self.keep_dominated = keep_dominated
        self.only_new_entries = only_new_entries
        self.hv_approx_samples = hv_approx_samples
        self.num_refpoints = num_refpoints


class MOArchive(Archive, ConfigurableObject):
    config_class = MOArchiveConfig

    def __init__(self, config: MOArchiveConfig, env: gym.Env, **kwargs) -> None:
        ConfigurableObject.__init__(self, config)
        Archive.__init__(self, max_size=config.max_size)
        self._archived_entries: DirtyList[MOArchiveEntry] = DirtyList([])
        self._hv: float = 0  # Empty archive has hypervolume of 0

        self.ideal_point, self.nadir_point, self.func = get_ideal_nadir(env)
        self.ref_dirs = get_ref_dirs(self.func.num_objectives, config.num_refpoints)
        self.num_refpoints: int = len(self.ref_dirs)

        self.secondary_criterion = self.secondary_criterion_class(
            ref_dirs=self.ref_dirs,
            ideal_point=self.ideal_point,
            nadir_point=self.nadir_point,
        )

    @property
    def archived_entries(self) -> DirtyList[MOArchiveEntry]:
        return self._archived_entries

    @archived_entries.setter
    def archived_entries(self, entries: list[MOArchiveEntry]) -> None:
        self._archived_entries = DirtyList(entries)
        self._hv = np.nan  # Reset the hypervolume score when the archive changes

    def _add(self, entries: list[ArchiveEntry]) -> list[dict[str, Any]]:
        """
        Add an entry to the archive. This method is meant to be implemented by subclasses.
        """
        self._hv = np.nan  # Reset the hypervolume score when the archive changes
        assert all(
            isinstance(entry, MOArchiveEntry) for entry in entries
        ), "All entries must be of type MOArchiveEntry"
        mo_entries = cast(list[MOArchiveEntry], entries)
        if self.only_new_entries:
            # Only new entries are added to the archive, so we replace the archived entries with the new ones
            self.archived_entries = mo_entries
        else:
            # All entries are added to the archive, so we extend the archived entries with the new ones
            self.archived_entries.extend(mo_entries)

        if self.keep_dominated != KeepDominated.ALL:
            # Only keep non-dominated entries in the archive
            self.archived_entries = MOArchive.remove_dominated(
                self.archived_entries, self.keep_dominated
            )

        if self.max_size is not None and len(self.archived_entries) > self.max_size:
            # We need to remove some entries from the archive to keep it within the max size
            max_rank = MOArchive.crit_rank(self.archived_entries, self.max_size)
            # Remember the critical entries as candidates
            critical_entries = [
                entry for entry in self.archived_entries if entry.rank == max_rank
            ]
            # For now, only add the safe entries to the archive, and remove the critical entries. We will add the critical entries back later based on the secondary criterion.
            self.archived_entries = [
                entry for entry in self.archived_entries if entry.rank < max_rank
            ]

            secondary_scores = self.secondary_criterion.score(
                crit_entries=critical_entries, accepted_entries=self.archived_entries
            )
            for entry, score in zip(critical_entries, secondary_scores):
                entry.secondary_score = score
            remaining_slots = self.max_size - len(self.archived_entries)
            last_entries = sorted(
                critical_entries, key=lambda e: e.secondary_score, reverse=True
            )[:remaining_slots]
            self.archived_entries.extend(last_entries)
        return [self.get_archive_stats()]

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
        non_dominated_mask = is_nondominated(
            points, keep_weakly=keep_weakly, maximise=False
        )
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
        for i, entry in enumerate(entries):
            entry.rank = int(ranks[i])
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
        if len(self.archived_entries) == 0:
            self._hv = 0.0
        if not np.isnan(self._hv) and not self.archived_entries.dirty:
            return self._hv
        points = np.array([entry.objectives for entry in self.archived_entries])
        if self.hv_approx_samples is None:
            hv = hypervolume(points, ref=self.nadir_point, maximise=False)
        else:
            hv = hv_approx(
                points,
                ref=self.nadir_point,
                maximise=False,
                nsamples=self.hv_approx_samples,
                seed=None,
                method="Rphi-FWE+",
            )
        self._hv = hv
        self.archived_entries.mark_clean()
        return hv

    def get_archive_stats(self) -> dict[str, Any]:
        """
        Return a dictionary with the current archive stats.
        """
        points = np.array([entry.objectives for entry in self.archived_entries])
        niches, dist_to_niches, _ = associate_to_niches(
            points, self.ref_dirs, self.ideal_point, self.nadir_point
        )
        niche_count = calc_niche_count(len(self.ref_dirs), niches)
        filled_niches = np.sum(niche_count > 0)
        for entry, niche, dist in zip(self.archived_entries, niches, dist_to_niches):
            entry.niche = int(niche)
            entry.dist_to_ref = float(dist)

        distance_to_ideal = np.linalg.norm(points - self.ideal_point, axis=1)
        filled_niche_distances = {}
        for niche in range(len(self.ref_dirs)):
            distance_to_ideal_niche = distance_to_ideal[niches == niche]
            avg_dist = (
                np.mean(distance_to_ideal_niche)
                if len(distance_to_ideal_niche) > 0
                else np.nan
            )
            min_dist = (
                np.min(distance_to_ideal_niche)
                if len(distance_to_ideal_niche) > 0
                else np.nan
            )
            sd_dist = (
                np.std(distance_to_ideal_niche)
                if len(distance_to_ideal_niche) > 0
                else np.nan
            )
            filled_niche_distances[niche] = {
                "avg": avg_dist,
                "min": min_dist,
                "std": sd_dist,
            }
        niche_perf_avg = np.mean(
            [
                filled_niche_distances[n]["avg"]
                for n in filled_niche_distances
                if not np.isnan(filled_niche_distances[n]["avg"])
            ]
        )
        niche_min_sd = np.std(
            [
                filled_niche_distances[n]["min"]
                for n in filled_niche_distances
                if not np.isnan(filled_niche_distances[n]["min"])
            ]
        )

        return {
            "size": self.size,
            "score": self.score,
            "avg_dist_to_niches": np.mean(dist_to_niches),
            "avg_niche_count": np.mean(niche_count),
            "std_niche_count": np.std(niche_count),
            "filled_niches": filled_niches,
            "avg_dist_to_ideal": np.mean(distance_to_ideal),
            "coverage": (
                filled_niches / len(self.ref_dirs) if len(self.ref_dirs) > 0 else 0.0
            ),
            "niche_perf_avg": niche_perf_avg,
            "niche_min_sd": niche_min_sd,
        }

    def get_all(self) -> list[ArchiveEntry]:
        """
        Return all entries in the archive.
        """
        all_entries: list[ArchiveEntry] = cast(
            list[ArchiveEntry], self.archived_entries
        )
        return all_entries

    def get(self, index: int) -> ArchiveEntry | None:
        """
        Return the entry at the given index in the archive.
        """
        if index < 0 or index >= len(self.archived_entries):
            return None
        return self.archived_entries[index]
