from matplotlib.artist import get
from jaix.env.utils.archive.rv_binning_strategy import (
    RVBinningStrategy,
    RVBinningStrategyConfig,
)
import numpy as np
from pymoo.core.individual import Individual


def test_get_bin():
    # Create a simple reference direction and ideal/nadir points
    ref_dirs = np.array([[0, 1], [1, 0]], dtype=float)
    ideal = np.array([0, 0], dtype=float)
    nadir = np.array([1, 1], dtype=float)

    config = RVBinningStrategyConfig()
    binning_strategy = RVBinningStrategy(config, ref_dirs, ideal, nadir)

    for _ in range(10):
        # Create a random Individual with objectives in the range [0, 1]
        sample = Individual()
        sample.F = np.random.rand(2) * 4  # Scale to [0, 4] for testing
        bin_index = binning_strategy.get_bin(sample)
        if sample.F[0] < sample.F[1]:
            assert bin_index == 0, f"Expected bin index 0, got {bin_index}"
        else:
            assert bin_index == 1, f"Expected bin index 1, got {bin_index}"


def test_get_k_nearest_bins():
    # Create a simple reference direction and ideal/nadir points
    ref_dirs = np.array(
        [[0, 1], [0.25, 0.75], [0.5, 0.5], [0.75, 0.25], [1, 0]], dtype=float
    )
    ideal = np.array([0, 0], dtype=float)
    nadir = np.array([1, 1], dtype=float)

    config = RVBinningStrategyConfig()
    binning_strategy = RVBinningStrategy(config, ref_dirs, ideal, nadir)

    # Test nearest bins for the first bin (index 0)
    nearest_bins = binning_strategy.get_k_nearest_bins(2, k=2)
    assert set(nearest_bins) == {
        1,
        3,
    }, f"Expected nearest bins [1, 3], got {nearest_bins}"
