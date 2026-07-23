from abc import ABC, abstractmethod
import pandas as pd
from typing import Any, Dict, List, Tuple
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes


class Archive(ABC):
    def __init__(self):
        # Record stats over time
        self.stats_rows = []  # List to store stats for each call
        self._stats = pd.DataFrame()

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

    @abstractmethod
    def get_archive_stats(self) -> Dict[str, Any]:
        """
        Return a dictionary with the current archive stats
        """
        pass

    @property
    @abstractmethod
    def score(self) -> float:
        """
        Return the score of the archive as a float
        """
        pass

    def simulate_add(self, sample: Any, fitness: float, **kwargs) -> Tuple[bool, float]:
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
    def _add(self, sample: Any, fitness: float) -> Dict[str, Any]:
        """
        Internal method to add a sample to the archive
        Returns a dictionary with the result of the addition
        """
        pass

    def add(self, sample: Any, fitness: float) -> float:
        """
        Add a sample to the archive and return the reward obtained from adding it
        """
        prev_score = self.score
        result_dict = self._add(sample, fitness)
        self.stats_rows.append(result_dict)
        new_score = self.score
        reward = new_score - prev_score
        return reward

    @abstractmethod
    def get_all(self) -> List[Any]:
        """
        Return all samples in the archive
        """
        pass

    def plot_stats(
        self,
        stat_names: List[str] | None = None,
        fig_path: str | None = None,
    ) -> Tuple[Figure, Axes]:
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
