import numpy as np
from pymoo.algorithms.moo.nsga3 import associate_to_niches
from sklearn.neighbors import KDTree
from ttex.config import Config, ConfigurableObject

from jaix.env.utils.archive.binning_strategy import BinningStrategy


class RVBinningStrategyConfig(Config):
    def __init__(self):
        super().__init__()


class RVBinningStrategy(BinningStrategy[np.ndarray], ConfigurableObject):
    config_class = RVBinningStrategyConfig

    def __init__(
        self,
        config: RVBinningStrategyConfig,
        ref_dirs: np.ndarray,
        ideal: np.ndarray,
        nadir: np.ndarray,
    ):
        ConfigurableObject.__init__(self, config)
        self.ref_dirs = ref_dirs
        self.ideal = ideal
        self.nadir = nadir
        self.kdtree = KDTree(
            ref_dirs
        )  # Build a KDTree for efficient nearest neighbor search

    def get_bin(self, sample: np.ndarray) -> int:
        """
        Given an objective vector, return the bin index it belongs to.
        """
        niches, _, _ = associate_to_niches(
            sample.reshape(1, -1), self.ref_dirs, self.ideal, self.nadir
        )
        return niches[0]  # Return the index of the niche/bin

    def get_k_nearest_bins(self, bidx: int, k: int) -> list:
        """
        Given a bin index, return the indices of the k nearest bins.
        """
        assert 0 <= bidx < len(self.ref_dirs), "Bin index out of range."
        # k must be less than or equal to the number of training points minus 1 (to exclude the point itself)
        k = min(k, len(self.ref_dirs) - 1)
        _, indices = self.kdtree.query(
            self.ref_dirs[bidx].reshape(1, -1), k=k + 1
        )  # +1 to include the bin itself
        nearest_bins = indices[0][1:]  # Exclude the bin itself
        return nearest_bins.tolist()
