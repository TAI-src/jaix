from ttex.config import Config, ConfigurableObject, ConfigurableObjectFactory as COF
import numpy as np
from typing import Any, Dict, List, Tuple, Type, Set
from jaix.env.utils.archive.archive import Archive
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
from jaix.env.utils.archive.binning_strategy import BinningStrategy


class BinArchiveConfig(Config):
    def __init__(
        self,
        n_bins: int,
        max_fitness: float,  # maximum fitness value for normalisation
        binning_strategy: Type[BinningStrategy],
        binning_config: Config,
        np_bin: int = 1,
        coverage_weight: float = 0.5,  # weight for coverage in the score function
    ):
        Config.__init__(self)
        self.n_bins = n_bins
        self.np_bin = np_bin
        self.coverage_weight = coverage_weight
        self.max_fitness = max_fitness
        self.binning_strategy = binning_strategy
        self.binning_config = binning_config

        assert np_bin >= 1, "np_bin must be at least 1"


class BinArchive(ConfigurableObject, Archive):
    config_class = BinArchiveConfig

    def __init__(self, config: BinArchiveConfig, **kwargs):
        Archive.__init__(self)
        ConfigurableObject.__init__(self, config)
        self.map: Dict[int, List[Tuple[float, float]]] = {
            bin_idx: [] for bin_idx in range(self.n_bins)
        }
        self.n_points = 0
        self.covered_bins = 0
        self.bin_set: Set[int] = set()  # Set to keep track of covered bins
        self.hit_counter = np.zeros(self.n_bins, dtype=int)  # Counter for hits per bin
        self.add_counter = np.zeros(
            self.n_bins, dtype=int
        )  # Counter for additions per bin
        self.replace_counter = np.zeros(
            self.n_bins, dtype=int
        )  # Counter for replacements per bin
        self.total_fitness = 0.0  # Total fitness of all samples in the archive
        self.binner = COF.create(self.binning_strategy, self.binning_config, **kwargs)

    @property
    def coverage(self) -> float:
        """
        Return the coverage of the archive as a float between 0 and 1
        """
        return self.covered_bins / self.n_bins

    @property
    def fitness(self) -> float:
        """
        Return the average fitness of the archive as a float
        """
        if self.n_points == 0:
            return np.nan
        return (
            self.total_fitness / self.n_points / self.max_fitness
        )  # Normalised fitness

    @property
    def score(self) -> float:
        """
        Return the score of the archive as a float
        Score is a weighted sum of coverage and average fitness
        """
        return (
            self.coverage_weight * self.coverage
            + (1 - self.coverage_weight) * self.fitness
        )

    def get_archive_stats(self, bin_stats: bool = False, **kwargs) -> Dict[str, Any]:
        """
        Return a dictionary with the current archive stats
        """
        stats_dict = {
            "n_points": self.n_points,
            "coverage": self.coverage,
            "fitness": self.fitness,
            "score": self.score,
            "covered_bins": self.covered_bins,
            "hit_counter": np.sum(self.hit_counter),
            "add_counter": np.sum(self.add_counter),
            "replace_counter": np.sum(self.replace_counter),
            "add_rate": np.sum(self.add_counter) / (np.sum(self.hit_counter) + 1e-8),
            "replace_rate": np.sum(self.replace_counter)
            / (np.sum(self.hit_counter) + 1e-8),
        }
        # Add the number of bins to the keys of the stats_dict
        stat_dict = {f"{key}_{self.n_bins}": value for key, value in stats_dict.items()}

        if bin_stats:
            # Add bin-specific stats to the stats rows
            pbin_fit = {
                f"fbin_{bidx}": (
                    np.min([fit for _, fit in self.map[bidx]])
                    if len(self.map[bidx]) > 0
                    else np.nan
                )
                for bidx in range(self.n_bins)
            }
            stat_dict.update(pbin_fit)

        for key, value in kwargs.items():
            stat_dict[key] = value

        return stat_dict

    def _add(self, sample: Any, fitness: float) -> Dict[str, Any]:
        """
        Add a sample to the archive if it is better than the current best in the bin.
        If the bin is empty, add the sample directly.
        If the bin is full, replace the worst sample if the new sample is better.
        Assuming lower fitness is better
        """
        bin_idx = self.binner.get_bin(sample)
        self.hit_counter[bin_idx] += 1  # Increment hit counter for this bin

        if len(self.map[bin_idx]) < self.np_bin:
            # bin has space, add the sample
            added = self._append(sample, fitness, bin_idx)
        else:
            added = self._replace(sample, fitness, bin_idx)
        # Update stats after adding/replacing
        stats = self.get_archive_stats(bin_stats=True, hit_bin=bin_idx, added=added)
        return stats

    def _append(self, sample: Any, fitness: float, bin_idx: int) -> bool:
        """
        Logic for adding a sample to bin that still has space
        """
        assert len(self.map[bin_idx]) < self.np_bin, "bin is full, cannot add sample"
        # bin is not full, add the sample
        self.map[bin_idx].append((sample, fitness))
        self.n_points += 1
        if bin_idx not in self.bin_set:
            self.covered_bins += 1
            self.bin_set.add(bin_idx)
        self.add_counter[bin_idx] += 1  # Increment add counter for this bin
        self.total_fitness += fitness  # Update total fitness
        return True

    def _replace(self, sample: Any, fitness: float, bin_idx: int) -> bool:
        """
        Logic for replacing the worst sample in a full bin
        """
        assert (
            len(self.map[bin_idx]) == self.np_bin
        ), "bin is not full, cannot replace sample"
        # bin is full, check if we can replace the worst sample
        worst_sample, worst_fitness = max(self.map[bin_idx], key=lambda x: x[1])
        if fitness < worst_fitness:
            # Replace the worst sample with the new one
            self.map[bin_idx].remove((worst_sample, worst_fitness))
            self.map[bin_idx].append((sample, fitness))
            self.replace_counter[bin_idx] += 1  # Increment replace counter for this bin
            # Update total fitness
            self.total_fitness += fitness - worst_fitness
            return True
        else:
            return False

    def get_elite(self, bin_idx: int) -> Tuple[Any, float]:
        """
        Return the best sample and its fitness in the given bin
        """
        if len(self.map[bin_idx]) == 0:
            return None, np.nan
        # Return the sample with the best fitness in the bin
        best_sample, best_fitness = min(self.map[bin_idx], key=lambda x: x[1])
        return best_sample, best_fitness

    def get_all(self) -> List[Tuple[Any, float]]:
        """
        Return all samples in the archive as a list of tuples (sample, fitness)
        """
        all_samples = []
        for bin_samples in self.map.values():
            all_samples.extend(bin_samples)
        return all_samples

    def plot_stats(
        self, stat_names: List[str] | None = None, fig_path: str | None = None
    ) -> Tuple[Figure, Axes]:
        """
        Plot the stats over time
        """
        if stat_names is None:
            stat_names = ["coverage", "add_rate", "replace_rate", "fitness", "score"]
            stat_names = [f"{col}_{self.n_bins}" for col in stat_names]
        return super().plot_stats(stat_names=stat_names, fig_path=fig_path)

    def plot_pbin_stats(self, fig_path: str | None = None) -> Tuple[Figure, Axes]:
        """
        Plot the number of hits per bin, split by successful additions
        """
        df = self.stats[["hit_bin", "added"]].copy()
        hits = df.groupby("hit_bin").size().reindex(range(self.n_bins), fill_value=0)
        adds = (
            df[df["added"]]
            .groupby("hit_bin")
            .size()
            .reindex(range(self.n_bins), fill_value=0)
        )
        fails = hits - adds
        # plot adds and fails as a stacked bar chart
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.bar(range(self.n_bins), adds, label="Successful adds", color="green")
        ax.bar(
            range(self.n_bins),
            fails,
            bottom=adds,
            label="Failed hits",
            color="red",
        )
        ax.set_xlabel("bin Index")
        ax.set_ylabel("Number of Hits")
        ax.set_title("bin Additions and failed hits")
        ax.legend()
        if fig_path is not None:
            plt.savefig(fig_path)
        return fig, ax

    def plot_pbin_history(self, fig_path: str | None = None) -> Tuple[Figure, Axes]:
        """
        Plot fitness and hits over time as a heatmap
        """
        bin_cols = [f"fbin_{bidx}" for bidx in range(self.n_bins)]
        n_iter = len(self.stats)
        fig, ax = plt.subplots(figsize=(len(bin_cols), n_iter))
        # Plot fitness heatmap
        sns.heatmap(
            self.stats[bin_cols],
            ax=ax,
            cmap="viridis",
            cbar_kws={"label": "Fitness"},
            xticklabels=False,
            yticklabels=False,
            linewidths=0,
            square=True,
        )
        # Add a box around the hit bin for each iteration
        # Colour depends on success
        df = self.stats[["hit_bin", "added"]].copy()
        for i, (bin_idx, added) in df.iterrows():
            color = "black" if added else "white"
            ax.add_patch(
                Rectangle(
                    (bin_idx, i),
                    1,
                    1,
                    fill=False,
                    edgecolor=color,
                    lw=2,
                )
            )
        ax.set_xlabel("bin Index")
        ax.set_ylabel("Iteration")
        ax.set_title("bin Fitness and Hits Over Time")
        if fig_path is not None:
            plt.savefig(fig_path)
        return fig, ax
