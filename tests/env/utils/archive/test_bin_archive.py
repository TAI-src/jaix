from jaix.env.utils.archive.bin_archive import BinArchive, BinArchiveConfig
from ttex.config import ConfigurableObjectFactory as COF, Config, ConfigurableObject
from jaix.env.utils.archive.binning_strategy import BinningStrategy
import numpy as np
import pytest
import os.path as osp
import random


class DummyBinningStrategyConfig(Config):
    def __init__(self, n_bins=5):
        Config.__init__(
            self,
        )
        self.n_bins = n_bins
        self.lower_bound = 0.0
        self.upper_bound = 10.0
        self.binning_strategy = "linear"


class DummyBinningStrategy(BinningStrategy, ConfigurableObject):
    config_class = DummyBinningStrategyConfig

    def __init__(self, config):
        ConfigurableObject.__init__(self, config)

    def get_bin(self, sample):
        # For testing purposes, we just return a fixed bin index
        return 0

    def get_k_nearest_bins(self, bidx, k):
        # For testing purposes, return k random bins
        bin_list = range(self.n_bins)
        sampled = random.sample(bin_list, k)
        # remove original bin index if present
        if bidx in sampled:
            sampled.remove(bidx)
        return sampled


def get_archive(pre_fill=False):
    config = BinArchiveConfig(
        n_bins=5,
        max_fitness=10.0,
        binning_strategy=DummyBinningStrategy,
        binning_config=DummyBinningStrategyConfig(),
        np_bin=2,
        coverage_weight=0.7,
    )
    archive = COF.create(BinArchive, config)
    if pre_fill:
        for i in range(5 * 2):  # Fill each bin with 2 samples
            sample = np.array([i, i + 1])
            fitness = float(10 - i)
            archive.add(sample, fitness)
    return archive


def test_init():
    archive = get_archive()
    assert archive.coverage == 0.0, "Initial coverage should be 0.0"
    assert np.isnan(archive.fitness), "Initial fitness should be NaN"
    assert np.isnan(archive.score), "Initial score should be NaN"
    archive_stats = archive.get_archive_stats()
    assert archive_stats["coverage_5"] == 0.0, "Stats coverage should be 0.0"


def test_append():
    archive = get_archive()
    sample = np.array([1.0, 2.0])
    fitness = 5.0
    added = archive._append(sample, fitness, bin_idx=0)
    assert added, "Sample should be added to the archive"
    x, f = archive.map[0][0]
    assert np.array_equal(x, sample), "Sample in archive should match the added sample"
    assert f == fitness, "Fitness in archive should match the added fitness"
    added = archive._append(sample, fitness, bin_idx=0)
    assert added, "Sample should be added to the archive even if it already exists"
    added = archive._append(sample, fitness, bin_idx=1)
    assert added, "Sample should be added to a different bin"
    assert len(archive.map[1]) == 1, "There should be one sample in bin 1"
    with pytest.raises(AssertionError):
        archive._append(sample, fitness, bin_idx=0)


def test_replace():
    archive = get_archive()
    sample1 = np.array([1.0, 2.0])
    fitness1 = 5.0
    archive._append(sample1, fitness1, bin_idx=0)
    sample2 = np.array([3.0, 4.0])
    fitness2 = 7.0
    with pytest.raises(AssertionError):
        archive._replace(sample2, fitness2, bin_idx=0)

    archive._append(sample1, fitness1, bin_idx=0)
    replaced = archive._replace(sample2, fitness2, bin_idx=0)
    assert not replaced, "Sample should not replace since it has worse fitness"
    assert (
        len(archive.map[0]) == archive.np_bin
    ), "There should be only 2 samples in bin 0 after replacement attempt"
    replaced = archive._replace(sample2, 3.0, bin_idx=0)
    assert replaced, "Sample should replace since it has better fitness"
    assert (
        archive.map[0][1][1] == 3.0
    ), "Fitness in archive should be updated to the new fitness"


def test_get_elite():
    archive = get_archive(pre_fill=True)
    _, fit = archive.get_elite(0)
    assert fit == 1.0, "Elite fitness in bin 0 should be 1.0"


def test_get_all():
    archive = get_archive(pre_fill=True)
    all_samples = archive.get_all()
    assert (
        len(all_samples) == 2
    ), "There should be 2 samples in the archive (only one bin gets filled)"


def test_plot_stats(tmp_path):
    archive = get_archive(pre_fill=True)
    plot_path = osp.join(tmp_path, "bin_stats.png")
    fig, ax = archive.plot_stats(fig_path=str(plot_path))
    assert fig is not None, "Figure should be created"
    assert ax is not None, "Axes should be created"
    assert osp.exists(plot_path), "Plot should be saved to the specified path"


def test_plot_pbin_stats(tmp_path):
    archive = get_archive(pre_fill=True)
    plot_path = osp.join(tmp_path, "pbin_stats.png")
    fig, ax = archive.plot_pbin_stats(fig_path=str(plot_path))
    assert fig is not None, "Figure should be created"
    assert ax is not None, "Axes should be created"
    assert osp.exists(plot_path), "Plot should be saved to the specified path"


def test_plot_pbin_history(tmp_path):
    archive = get_archive(pre_fill=True)
    plot_path = osp.join(tmp_path, "pbin_history.png")
    fig, ax = archive.plot_pbin_history(fig_path=str(plot_path))
    assert fig is not None, "Figure should be created"
    assert ax is not None, "Axes should be created"
    assert osp.exists(plot_path), "Plot should be saved to the specified path"


def test_get_closest_bin():
    archive = get_archive(pre_fill=True)
    sample = np.array([2.0, 3.0])
    bin = archive.binner.get_bin(sample)
    closest_bin = archive.get_closest_non_empty_bin(bin)
    assert closest_bin == bin, "Closest bin index should be 0 for the given sample"

    closest_bin = archive.get_closest_non_empty_bin(bin + 2)
    assert closest_bin == bin, "Closest bin index should be 0 for the given sample"


def test_closest_elite():
    archive = get_archive(pre_fill=True)
    sample = np.array([2.0, 3.0])
    bin = archive.binner.get_bin(sample)
    elite_sample, elite_fitness, closest_bin = archive.get_closest_elite(bin)
    assert closest_bin == bin, "Closest bin index should be 0 for the given sample"

    elite_sample, elite_fitness, closest_bin = archive.get_closest_elite(bin + 2)
    assert closest_bin == bin, "Closest bin index should be 0 for the given sample"
    assert (
        elite_fitness == archive.get_elite(bin)[1]
    ), "Elite fitness should match the fitness of the closest bin"
