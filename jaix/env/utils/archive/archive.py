from abc import ABC, abstractmethod

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from typing import TypeVar, Generic, Any

T = TypeVar("T")


class ArchiveEntry(ABC, Generic[T]):
    """
    A general archive entry class that stores all information about an entry in the archive.
    """

    @abstractmethod
    def parse(self) -> T:
        """
        Return a dictionary with the parsed information from the archive entry
        This needs to be implemented by the subclass, and should return a dictionary with the information needed for the archive to store the entry.
        """


class Archive(ABC):
    """
    A general archive class that can be used to store samples and their fitness values.
    This mainly implements recording and plotting stats over time, but the actual archive implementation is left to the subclasses.
    """

    def __init__(self, max_size: int | None = None):
        # Record stats over time
        self._max_size = max_size
        self.reset()

    def reset(self) -> None:
        self.stats_rows: list[dict] = []
        self._stats: pd.DataFrame = pd.DataFrame()

    @property
    def stats(self) -> pd.DataFrame:
        """
        Update the stats dataframe with the current stats rows
        and return
        """
        if len(self.stats_rows) == 0:
            return pd.DataFrame()
        if len(self.stats_rows) != len(self._stats):
            self._stats = pd.DataFrame(self.stats_rows)
        return self._stats

    @property
    def max_size(self) -> int | None:
        return self._max_size

    @property
    @abstractmethod
    def size(self) -> int:
        """
        Return the current size of the archive
        """

    @abstractmethod
    def get_archive_stats(self) -> dict[str, Any]:
        """
        Return a dictionary with the current archive stats
        """

    @property
    @abstractmethod
    def score(self) -> float:
        """
        Return the score of the archive as a float
        """

    def simulate_add(self, sample: Any, fitness: float, **kwargs) -> tuple[bool, float]:
        """
        Simulate adding a sample to the archive without actually adding it
        Returns a tuple (added, reward) where added is a boolean
        indicating if the sample would be added and
        reward is the reward that would be obtained from adding the sample
        """
        raise NotImplementedError(
            "simulate_add is not implemented for this archive type"
        )

    @abstractmethod
    def _add(self, entry: ArchiveEntry) -> dict[str, Any]:
        """
        Internal method to add a sample to the archive
        Returns a dictionary with the result of the addition
        """

    def add(self, entry: ArchiveEntry) -> float:
        """
        Add a sample to the archive and return the reward obtained from adding it
        """
        prev_score = self.score
        result_dict = self._add(entry)
        self.stats_rows.append(result_dict)
        new_score = self.score
        reward = new_score - prev_score
        return reward

    @abstractmethod
    def get_all(self) -> list[ArchiveEntry]:
        """
        Return all samples in the archive
        """

    @abstractmethod
    def get(self, index: int) -> ArchiveEntry | None:
        """
        Return the archive entry at the given index as (sample, fitness), or None if unavailable.
        """

    def plot_stats(
        self,
        stat_names: list[str] | None = None,
        fig_path: str | None = None,
    ) -> tuple[Figure, Axes]:
        """
        Plot the stats over time
        """

        if stat_names is None:
            stat_names = self.stats.columns.tolist()
        fig, ax = plt.subplots(figsize=(10, 6))
        self.stats[stat_names].plot(ax=ax)
        ax.set_title("Archive Stats Over Time")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Value")
        ax.grid()
        ax.legend(title="Stats")
        if fig_path is not None:
            plt.savefig(fig_path)
        return fig, ax
