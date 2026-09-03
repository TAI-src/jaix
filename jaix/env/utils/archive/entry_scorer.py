from abc import ABC, abstractmethod

import numpy as np
from moocore import hv_contributions
from pymoo.algorithms.moo.nsga3 import associate_to_niches, calc_niche_count, niching

from jaix.env.utils.archive.archive import ArchiveEntry


class EntryScorer(ABC):

    @abstractmethod
    def score(
        self, crit_entries: list[ArchiveEntry], accepted_entries: list[ArchiveEntry]
    ) -> list[float]:
        """
        Score the critical individuals according to the secondary criterion.
        The higher the score, the better the individual is considered.
        """
        ...


class HVContributionScorer(EntryScorer):
    def __init__(self, nadir_point: np.ndarray, **kwargs) -> None:
        self.nadir_point = nadir_point

    def score(
        self, crit_entries: list[ArchiveEntry], accepted_entries: list[ArchiveEntry]
    ) -> list[float]:
        points = np.array(
            [getattr(entry, "objectives", None) for entry in crit_entries]
        )
        scores = hv_contributions(
            points, ref=self.nadir_point, maximise=False, ignore_dominated=True
        )
        return scores.tolist()


class ReferenceVectorDistanceScorer(EntryScorer):
    def __init__(
        self,
        ref_dirs: np.ndarray,
        ideal_point: np.ndarray,
        nadir_point: np.ndarray,
        **kwargs,
    ) -> None:
        self.ref_dirs = ref_dirs
        self.ideal_point = ideal_point
        self.nadir_point = nadir_point

    def score(
        self, crit_entries: list[ArchiveEntry], accepted_entries: list[ArchiveEntry]
    ) -> list[float]:
        c_points = np.array(
            [getattr(entry, "objectives", None) for entry in crit_entries]
        )
        c_niches, c_dist_to_niches, _ = associate_to_niches(
            c_points, self.ref_dirs, self.ideal_point, self.nadir_point
        )
        if len(accepted_entries) > 0:
            accepted_points = np.array(
                [getattr(entry, "objectives", None) for entry in accepted_entries]
            )
            accepted_niches, _, _ = associate_to_niches(
                accepted_points, self.ref_dirs, self.ideal_point, self.nadir_point
            )
            niche_count = calc_niche_count(len(self.ref_dirs), accepted_niches)
        else:
            niche_count = np.zeros(len(self.ref_dirs), dtype=int)

        indices = niching(
            list(range(len(crit_entries))),
            len(crit_entries),
            niche_count,
            c_niches,
            c_dist_to_niches,
        )
        # The score of the entry is the inverse position in the sorted list of indices, so that the first entry has the highest score
        scores = [len(indices) - indices.index(i) for i in range(len(crit_entries))]
        return scores
