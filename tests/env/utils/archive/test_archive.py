import os.path as osp
from typing import Any

import gymnasium as gym
from jaix.env.utils.archive.archive import Archive, ArchiveEntry


class DummyArchiveEntry(ArchiveEntry[dict[str, Any]]):
    def __init__(self, sample: Any, fitness: float):
        self.sample = sample
        self.fitness = fitness

    def parse(self) -> dict[str, Any]:
        return {"sample": self.sample, "fitness": self.fitness}


class DummyArchive(Archive):
    def __init__(self, max_size: int, **kwargs):
        super().__init__(max_size=max_size)
        self._score = 0.0
        self.num_points = 0
        self.added_samples = []

    @property
    def archive_entry_type(self) -> type[ArchiveEntry]:
        return DummyArchiveEntry

    @property
    def score(self) -> float:
        return self._score

    def get_archive_stats(self) -> dict[str, Any]:
        return {"score": self._score, "num_points": self.num_points}

    def _add(self, entry: ArchiveEntry) -> dict[str, Any]:
        assert isinstance(entry, DummyArchiveEntry), "Entry must be a DummyArchiveEntry"
        # Simple archive of max_size, we just remove the oldest sample if we exceed max_size
        if self.max_size is not None and self.num_points >= self.max_size:
            self.added_samples.pop(0)
        self.added_samples.append(entry)
        self.num_points += 1
        self._score += entry.fitness
        return self.get_archive_stats()

    @property
    def size(self) -> int:
        return len(self.added_samples)

    def get_all(self):
        return self.added_samples

    def get(self, index: int):
        if index < 0 or index >= len(self.added_samples):
            return None
        else:
            return self.added_samples[index]


def test_archive_add():
    archive = DummyArchive(max_size=10)
    initial_score = archive.score
    initial_num_points = archive.num_points

    # Add a sample with fitness 10.0
    reward = archive.add(DummyArchiveEntry(sample="sample1", fitness=10.0))
    assert reward == 10.0, "Reward should be equal to the fitness value"
    assert archive.score == initial_score + 10.0, "Score should be updated correctly"
    assert (
        archive.num_points == initial_num_points + 1
    ), "Number of points should be incremented"

    # Add another sample with fitness 5.0
    reward = archive.add(DummyArchiveEntry(sample="sample2", fitness=5.0))
    assert reward == 5.0, "Reward should be equal to the fitness value"
    assert archive.score == initial_score + 15.0, "Score should be updated correctly"
    assert (
        archive.num_points == initial_num_points + 2
    ), "Number of points should be incremented"

    # Test the stats property
    stats = archive.stats
    assert (
        stats["score"].iloc[-1] == archive.score
    ), "Stats should reflect the current score"
    assert (
        stats["num_points"].iloc[-1] == archive.num_points
    ), "Stats should reflect the current number of points"


def test_plot(tmp_path):
    archive = DummyArchive(max_size=10)
    # Add some samples to the archive
    for i in range(5):
        archive.add(DummyArchiveEntry(sample=f"sample{i}", fitness=float(i) * 2))

    fig_path = osp.join(tmp_path, "archive_stats.png")
    fig, ax = archive.plot_stats(
        fig_path=fig_path
    )  # This should generate a plot without errors
    # Check if the figure is saved correctly
    assert osp.exists(fig_path), "Figure should be saved to the specified path"
    # Check if the figure and axes are returned correctly
    assert fig is not None, "Figure should not be None"
    assert ax is not None, "Axes should not be None"
