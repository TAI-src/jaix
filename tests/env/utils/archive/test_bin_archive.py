import os.path as osp
import random
from typing import Any

import numpy as np
import pytest
from jaix.env.utils.archive.bin_archive import (
    BinArchive,
    BinArchiveConfig,
    BinArchiveEntry,
)
from jaix.env.utils.archive.binning_strategy import BinningStrategy
from ttex.config import Config, ConfigurableObject
from ttex.config import ConfigurableObjectFactory as COF


class DummyArchiveEntry(BinArchiveEntry[Any]):
    def __init__(self, sample: Any, fitness: float):
        self.sample = sample
        self._fitness = fitness

    def parse(self) -> tuple[Any, float]:
        return self.sample, self._fitness


class DummyBinningStrategyConfig(Config):
    def __init__(self, n_bins=5):
        Config.__init__(
            self,
        )
        self.n_bins = n_bins
        self.lower_bound = 0.0
        self.upper_bound = 10.0
        self.binning_strategy = "linear"


class DummyBinningStrategy(BinningStrategy[Any], ConfigurableObject):
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


def get_archive(pre_fill=False, allow_close_elites=True):
    config = BinArchiveConfig(
        n_bins=5,
        max_fitness=10.0,
        binning_strategy=DummyBinningStrategy,
        binning_config=DummyBinningStrategyConfig(),
        archive_entry_class=DummyArchiveEntry,
        np_bin=2,
        coverage_weight=0.7,
        allow_close_elites=allow_close_elites,
    )
    archive = COF.create(BinArchive, config)
    if pre_fill:
        for i in range(5 * 2):  # Fill each bin with 2 samples
            sample = np.array([i, i + 1])
            fitness = float(10 - i)
            entry = DummyArchiveEntry(sample, fitness)
            a_stats = archive.get_archive_stats()
            prev_counter = a_stats["add_counter_5"] + a_stats["replace_counter_5"]
            archive.add(entry)
            assert (
                entry.bin_idx is not None
            ), "Entry should have a bin index after being added"
            a_stats = archive.get_archive_stats()
            new_counter = a_stats["add_counter_5"] + a_stats["replace_counter_5"]
            if new_counter > prev_counter:
                assert (
                    entry.added is True
                ), "Entry should be marked as added if it increased the archive size"
            else:
                assert (
                    entry.added is False
                ), "Entry should not be marked as added if it did not increase the archive size"
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
    entry = DummyArchiveEntry(sample=np.array([1.0, 2.0]), fitness=5.0)
    added = archive._append(entry, bin_idx=0)
    assert added, "Sample should be added to the archive"
    map_entry = archive.map[0][0]
    assert np.array_equal(
        map_entry.bin_sample, entry.bin_sample
    ), "Sample in archive should match the added sample"
    assert (
        entry.fitness == map_entry.fitness
    ), "Fitness in archive should match the added fitness"
    added = archive._append(entry, bin_idx=0)
    assert added, "Sample should be added to the archive even if it already exists"
    added = archive._append(entry, bin_idx=1)
    assert added, "Sample should be added to a different bin"
    assert len(archive.map[1]) == 1, "There should be one sample in bin 1"
    with pytest.raises(AssertionError):
        archive._append(entry, bin_idx=0)


def test_replace():
    archive = get_archive()
    entry_1 = DummyArchiveEntry(sample=np.array([1.0, 2.0]), fitness=5.0)
    archive._append(entry_1, bin_idx=0)
    entry_2 = DummyArchiveEntry(sample=np.array([3.0, 4.0]), fitness=7.0)
    with pytest.raises(AssertionError):
        archive._replace(entry_2, bin_idx=0)

    archive._append(entry_1, bin_idx=0)
    replaced = archive._replace(entry_2, bin_idx=0)
    assert not replaced, "Sample should not replace since it has worse fitness"
    assert (
        len(archive.map[0]) == archive.np_bin
    ), "There should be only 2 samples in bin 0 after replacement attempt"
    entry_3 = DummyArchiveEntry(sample=entry_2.sample, fitness=3.0)
    replaced = archive._replace(entry_3, bin_idx=0)
    assert replaced, "Sample should replace since it has better fitness"
    entry_4 = archive.map[0][1]
    assert (
        entry_4.fitness == 3.0
    ), "Fitness in archive should be updated to the new fitness"


def test_get_elite():
    archive = get_archive(pre_fill=True)
    entry = archive.get_elite(0)
    assert entry.fitness == 1.0, "Elite fitness in bin 0 should be 1.0"


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
    elite_entry, closest_bin = archive.get_closest_elite(bin)
    assert closest_bin == bin, "Closest bin index should be 0 for the given sample"

    elite_entry, closest_bin = archive.get_closest_elite(bin + 2)
    assert closest_bin == bin, "Closest bin index should be 0 for the given sample"
    assert (
        elite_entry.fitness == archive.get_elite(bin).fitness
    ), "Elite fitness should match the fitness of the closest bin"


def test_get():
    archive = get_archive(pre_fill=True)
    entry = archive.get(0)
    assert entry is not None, "Entry should not be None for valid index"

    entry = archive.get(2)
    assert (
        entry is not None
    ), "Entry should not be None for valid index if closest is allowed"

    entry = archive.get(100)
    assert entry is None, "entry should be None for out-of-range bin index"

    archive = get_archive(pre_fill=True, allow_close_elites=False)
    entry = archive.get(2)
    assert (
        entry is None
    ), "entry should be None for empty bin index when closest is not allowed"


def test_size():
    archive = get_archive(pre_fill=True)
    assert archive.size == archive.covered_bins
