from pymoo.algorithms.moo.nsga3 import associate_to_niches
from pymoo.core.individual import Individual
import numpy as np
from sklearn.neighbors import KDTree
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns
from ttex.config import ConfigurableObject, Config
from jaix.env.utils.decomp_archive.rec_vec_defaults import get_ref_dirs, get_pf

import logging
from jaix.utils.globals import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)


class RecombArchiveConfig(Config):
    def __init__(
        self,
        num_ref_points: int | str,
        n_per_bucket: int = 1,
        coverage_weight: float = 0.5,
    ):
        self.num_ref_points = num_ref_points
        self.n_per_bucket = n_per_bucket
        self.coverage_weight = coverage_weight

        assert self.n_per_bucket >= 1, "n_per_bucket must be at least 1"
        if isinstance(self.num_ref_points, int):
            assert self.num_ref_points >= 1, "num_ref_points must be at least 1"
        else:
            assert self.num_ref_points in [
                "original",
                "energy_small",
                "energy_medium",
                "energy_large",
            ], "num_ref_points must be one of 'original', 'energy_small', 'energy_medium', 'energy_large'"


class RecombArchive(ConfigurableObject):
    config_class = RecombArchiveConfig
    """
    This is the archive that keeps track of the best solution found for each bucket,
    i.e. niche determined by reference direction
    """

    def __init__(self, config):
        super().__init__(config)

        self.num_points = 0

        self.covered_niches = 0  # how many niches have at least one solution in them
        self.total_ideal_distance: float = float(
            0  # average distance to ideal point across all elites in the archive
        )
        self.niche_set = set()  # Set to keep track of covered niches
        # Record distance to ideal point for each niche for each iteration
        # columns are niches plus stats, rows are iterations
        self.stats_rows = []
        self._stats = pd.DataFrame()  # DataFrame to store stats over time

    def initialize(self, problem):
        """
        Initialize the archive based on the problem
        """
        pf, self.ideal, self.nadir = get_pf(
            problem
        )  # Using default number of refpoints for estimate
        self.max_improvement: float = float(
            np.linalg.norm(self.nadir - self.ideal)
        )  # maximum possible improvement (distance from nadir to ideal)
        assert self.max_improvement > 0, "Ideal and nadir points must be different"
        self.ref_dirs = get_ref_dirs(problem.n_obj, self.num_ref_points)
        self.kdtree = KDTree(
            self.ref_dirs
        )  # KDTree for efficient nearest neighbor search
        self.qd_scores = np.zeros(
            len(self.ideal)
        )  # Initialize QD scores for each objective
        self.hit_counter = np.zeros(
            len(self.ref_dirs)
        )  # Counter for the number of times each niche has been hit
        self.add_counter = np.zeros(
            len(self.ref_dirs)
        )  # Counter for the number of times each niche has been added to
        self.replace_counter = np.zeros(
            len(self.ref_dirs)
        )  # Counter for the number of times each niche has been replaced

        self.map = {niche_idx: [] for niche_idx in range(len(self.ref_dirs))}

    @property
    def stats(self) -> pd.DataFrame:
        """
        Update the stats dataframe with the current stats rows
        and return
        """
        if len(self.stats_rows) != len(self._stats):
            self._stats = pd.DataFrame(self.stats_rows)
        return self._stats

    def add(self, sol: Individual | np.ndarray, F: np.ndarray | None = None) -> float:
        if isinstance(sol, np.ndarray):
            assert F is not None, "If sol is a numpy array, F must be provided"
            sol = Individual(X=sol, F=F)
        elif isinstance(sol, Individual):
            if F is not None:
                sol.F = F  # Update the objective values if provided
            else:
                F = sol.F  # Use the objective values from the Individual
        else:
            raise ValueError("sol must be an Individual or a numpy array")
        assert F is not None, "Objective values F must be provided"

        niches, _, _ = associate_to_niches(
            F.reshape(1, -1), self.ref_dirs, self.ideal, self.nadir
        )
        niche = niches[0]  # Get the first (and only) niche index
        self.hit_counter[niche] += 1  # Increment the hit counter for the niche
        # Get distance to ideal point for sorting later
        dist_ideal = float(np.linalg.norm(F - self.ideal))
        if len(self.map[niche]) < self.n_per_bucket:
            # If the niche is not full, add the solution to the niche
            logger.debug(f"Adding solution to niche {niche} (not full)")
            added = self._add_sol_free(sol, F, dist_ideal, niche)
        else:
            # If the niche is full, check if new solution is replacing and do the replacement if so
            added = self._add_sol_occupied(sol, F, dist_ideal, niche)
        stats = self.get_archive_stats()
        stat_plots_dict = {
            f"niche_{niche_idx}": (
                self.map[niche_idx][0][2] if len(self.map[niche_idx]) > 0 else np.nan
            )
            for niche_idx in self.map
        }
        stat_plots_dict.update(stats)
        stat_plots_dict["hit_niche"] = niche
        stat_plots_dict["added"] = added
        self.stats_rows.append(stat_plots_dict)

        logger.debug(f"Archive statistics: {stats}")
        return self.score_archive(self.coverage_weight)

    def _add_sol_free(self, sol, F, dist_ideal, niche):
        """
        Add solution to the archive when the niche is not full.
        Update trackers for score calculation
        """
        assert len(self.map[niche]) < self.n_per_bucket, "Niche is already full"
        self.map[niche].append((sol, F, dist_ideal))
        self.num_points += 1  # Increment the number of points in the archive
        self.add_counter[niche] += 1  # Increment the add counter for the niche
        self.total_ideal_distance += dist_ideal  # Add the new distance to the total
        if len(self.map[niche]) == 1:
            self.covered_niches += (
                1  # Increment covered niches if this is the first solution in the niche
            )
        self.niche_set.add(niche)  # Add the niche to the set of covered niches
        assert (
            len(self.niche_set) == self.covered_niches
        ), "Mismatch in covered niches count"
        # sort the points in the niche by distance to ideal point
        self.map[niche].sort(key=lambda x: x[2])
        # Update QD scores for each objective
        # We are assuming that the objectives are to be minimized, so we take the negative of the fitness values
        self.qd_scores += -F  # Update QD scores for each objective
        return True

    def _add_sol_occupied(self, sol, F, dist_ideal, niche):
        """
        Add solution to the archive when the niche is full.
        Check if the new solution is better than the worst solution in the niche and replace if so
        Update trackers accordingly
        """
        # check the last point in the list (furthest from ideal point) and replace if new point is better
        _, old_F, old_dist = self.map[niche][-1]
        if dist_ideal < old_dist:
            logger.debug(
                f"Replacing solution in niche {niche} (full) with improvement {old_dist - dist_ideal}"
            )
            self.replace_counter[
                niche
            ] += 1  # Increment the replace counter for the niche
            self.map[niche][-1] = (sol, F, dist_ideal)
            self.total_ideal_distance -= (
                old_dist  # Remove the old distance from the total
            )
            self.total_ideal_distance += dist_ideal  # Add the new distance to the total
            # sort the points in the niche by distance to ideal point
            self.map[niche].sort(key=lambda x: x[2])
            # Update QD scores for each objective
            self.qd_scores += old_F  # Remove the old fitness values from the QD scores
            self.qd_scores += -F  # Add the new fitness values to the QD scores
            return True
        else:
            logger.debug(
                f"Discarding solution for niche {niche} (full) with no improvement {old_dist - dist_ideal}"
            )
            return False

    def score_archive(self, coverage_weight: float) -> float:
        """
        Score the archive based on coverage and average distance to ideal point
        """
        coverage = self.coverage()
        avg_dist_ideal = self.average_ideal_distance()
        score = (1 - coverage) * coverage_weight + avg_dist_ideal * (
            1 - coverage_weight
        )
        return score

    def coverage(self) -> float:
        """
        Get the coverage of the archive, i.e. the number of niches that have at least one solution
        """
        return self.covered_niches / len(self.map)

    def average_ideal_distance(self) -> float:
        """
        Get the average distance to the ideal point across all elites in the archive
        """
        if self.num_points == 0:
            return 1.0  # If there are no points, return maximum distance
        else:
            return self.total_ideal_distance / (self.num_points * self.max_improvement)

    def qd_score(self) -> np.ndarray:
        """
        Get the QD scores for the archive,
        i.e. one score per objective which is the sum of the fitness values (-1 for min)
        normalised across the number of points
        """
        if self.num_points == 0:
            return np.zeros(
                len(self.ideal)
            )  # If there are no points, return zero scores
        else:
            return self.qd_scores / self.num_points  # Normalize by the number of poins

    def get_archive_stats(self) -> dict:
        """
        Get the stats of the archive, including coverage, average distance to ideal point, and QD scores
        """
        stats = {
            f"coverage_{len(self.ref_dirs)}": self.coverage(),
            f"avg_dist_ideal_{len(self.ref_dirs)}": self.average_ideal_distance(),
            f"num_points_{len(self.ref_dirs)}": self.num_points,
            "hit_counter": np.sum(
                self.hit_counter
            ),  # Include the hit counter in the stats
            f"add_rate_{len(self.ref_dirs)}": np.sum(self.add_counter)
            / (np.sum(self.hit_counter) + 1e-8),
            f"replace_rate_{len(self.ref_dirs)}": np.sum(self.replace_counter)
            / (np.sum(self.hit_counter) + 1e-8),
        }
        stats.update(
            {
                f"qd_score_{len(self.ref_dirs)}_{i}": score
                for i, score in enumerate(self.qd_score())
            }
        )
        return stats

    def get_elite(self, niche: int) -> tuple[Individual, np.ndarray, float] | None:
        """Get the (sol, F, dist_ideal) tuple for a niche, or None if the niche is empty."""
        return self.map[niche][0] if len(self.map[niche]) > 0 else None

    def get_closest_elite(
        self,
        niche: int,  # Niche index to find the closest elite for
        safety_k: int = 5,  # Safety factor for the number of closest niches to consider
        max_k: int = -1,  # Maximum number of closest niches to consider (if -1, consider all)
    ) -> tuple[Individual, np.ndarray, float, int] | None:
        """
        Get the elite solution from the niche. If it doesn't exist, find the closest niche and return its elite solution.
        The number of niches to check is determined by the safety factor and the number of covered niches, but can be limited by max_k.
        If the neighboring niches don't have elites, return None.
        """
        elite = self.get_elite(niche)
        if elite is not None:
            return (*elite, niche)
        if self.covered_niches == 0:
            return None  # No elites in the archive
        if max_k == -1:
            max_k = len(self.ref_dirs)
        # Find the closest niche with an elite solution
        k = int(
            min(max(1, safety_k) * len(self.ref_dirs) / self.covered_niches, max_k)
        )  # Determine the number of closest niches to consider based on the safety factor and the number of covered niches
        logger.debug(
            f"Looking for closest elite to niche {niche} with k={k} given safety_k={safety_k}, max_k={max_k} and {self.covered_niches} covered niches"
        )
        _, indices = self.kdtree.query(
            self.ref_dirs[niche].reshape(1, -1), k=int(k)
        )  # Query the KDTree for the closest niches
        # remove indices based on the niche set to only consider niches that have elites
        indices = [idx for idx in indices[0] if idx in self.niche_set]
        if len(indices) == 0:
            return None  # No niches with elites found
        closest_niche = indices[0]  # Get the closest niche with an elite solution
        distance = np.linalg.norm(self.ref_dirs[niche] - self.ref_dirs[closest_niche])
        logger.debug(
            f"Closest niche to {niche} is {closest_niche} with distance {distance:.4f} found {len(indices)} niches with elites for k {k}"
        )
        closest_elite = self.get_elite(closest_niche)
        assert closest_elite is not None, "Closest niche should have an elite solution"
        return (
            *closest_elite,
            closest_niche,
        )  # Return the elite solution from the closest niche

    def get_point(self, niche: int) -> tuple[np.ndarray, int]:
        """
        Get the best point for a given niche, or None if the niche is empty
        """
        # TODO: should make this nicer at some point
        elite = self.get_closest_elite(niche)
        x = None
        if elite is None:
            # get a random point from the archive
            random_niche = np.random.choice(list(self.niche_set))
            distance = np.linalg.norm(
                self.ref_dirs[niche] - self.ref_dirs[random_niche]
            )
            logger.debug(
                f"No elite found for niche {niche}, returning random elite from niche {random_niche} with distance {distance:.4f}"
            )
            elite = self.get_elite(random_niche)
            assert elite is not None, "Random niche should have an elite solution"
            individual, _, _ = elite  # type: ignore
            x = individual.X
            ret_niche = random_niche
        else:
            individual, _, _, ret_niche = elite  # type: ignore
            x = individual.X
        assert x is not None, "Best point for the niche should not be None"
        assert isinstance(x, np.ndarray), "Returned point should be a numpy array"
        return (
            x,
            ret_niche,
        )  # Return the best point for the niche and the niche index of the elite solution

    def get_pop(self) -> list[Individual]:
        """
        Get the current population of elites in the archive, i.e. the best solution for each niche
        """
        return [self.map[niche][0][0] for niche in self.map if len(self.map[niche]) > 0]

    def plot_heatmap(self, fig_path: str | None = None):
        """
        Plot how the niches are filling up over time
        """

        # Plot the distance to ideal point for each niche over iterations in a heatmap
        niche_cols = [col for col in self.stats.columns if col.startswith("niche_")]
        iterations = len(self.stats)
        fig, ax = plt.subplots(figsize=(len(niche_cols), iterations))
        sns.heatmap(
            self.stats[niche_cols],
            ax=ax,
            cmap="viridis",
            cbar_kws={"label": "Distance to Ideal Point"},
            xticklabels=False,
            yticklabels=False,
            linewidths=0,
            square=True,
        )
        # Add a box around the cell that was indicated by hit_niche for each iteration
        # if the solution was added, add a black box, if it was not added, add a white box
        subset = self.stats[["hit_niche", "added"]]
        for iteration, (hit_niche, added) in subset.iterrows():
            color = "black" if added else "white"
            rect = Rectangle(
                (hit_niche, iteration),
                1,
                1,
                fill=False,
                edgecolor=color,
                linewidth=2,
            )
            ax.add_patch(rect)

        ax.set_title("Distance to Ideal Point for Each Niche Over Iterations")
        ax.set_ylabel("Iterations")
        ax.set_xlabel("Niches")
        if fig_path is not None:
            # Save the figure
            fig.savefig(fig_path, bbox_inches="tight")
        return fig, ax

    def plot_hits_add_by_niche(self, fig_path: str | None = None):
        """
        Plot the number of hits and adds for each niche
        """

        subset = self.stats[["hit_niche", "added"]]
        hits = subset.groupby("hit_niche").size()
        adds = subset[subset["added"]].groupby("hit_niche").size()
        unsuccessful_hits = hits - adds
        fig, ax = plt.subplots(figsize=(10, 6))
        # plot adds and unsuccessful hits as stacked bar chart
        ax.bar(
            hits.index,
            adds,
            label="Adds",
            color="blue",
        )
        ax.bar(
            hits.index,
            unsuccessful_hits,
            bottom=adds,
            label="Unsuccessful Hits",
            color="orange",
        )

        ax.set_title("Number of Hits and Adds for Each Niche")
        ax.set_ylabel("Count")
        ax.set_xlabel("Niche Index")
        ax.legend()
        if fig_path is not None:
            # Save the figure
            fig.savefig(fig_path, bbox_inches="tight")
        return fig, ax

    def plot_stats_time(self, fig_path: str | None = None):
        """
        Plot the coverage, add rate, and replace rate over time
        """
        stat_cols = [
            f"coverage_{len(self.ref_dirs)}",
            f"add_rate_{len(self.ref_dirs)}",
            f"replace_rate_{len(self.ref_dirs)}",
        ]
        fig, ax = plt.subplots(figsize=(10, 6))
        self.stats[stat_cols].plot(ax=ax)
        ax.set_title("Archive Statistics Over Time")
        ax.set_ylabel("Value")
        ax.set_xlabel("Iterations")
        ax.legend(title="Statistics")
        if fig_path is not None:
            # Save the figure
            fig.savefig(fig_path, bbox_inches="tight")
        return fig, ax
