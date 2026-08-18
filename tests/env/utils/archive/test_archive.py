from jaix.env.utils.archive.archive import Archive
from typing import Any, Dict
import os.path as osp


class DummyArchive(Archive):
    def __init__(self):
        super().__init__()
        self._score = 0.0
        self.num_points = 0

    @property
    def score(self) -> float:
        return self._score

    def get_archive_stats(self) -> Dict[str, Any]:
        return {"score": self._score, "num_points": self.num_points}

    def _add(self, sample: Any, fitness: float, **kwargs) -> Dict[str, Any]:
        # For testing purposes, we just increment the score by the fitness value
        self.num_points += 1
        self._score += fitness
        return self.get_archive_stats()

    def get_all(self):
        return []


def test_archive_add():
    archive = DummyArchive()
    initial_score = archive.score
    initial_num_points = archive.num_points

    # Add a sample with fitness 10.0
    reward = archive.add(sample="sample1", fitness=10.0)
    assert reward == 10.0, "Reward should be equal to the fitness value"
    assert archive.score == initial_score + 10.0, "Score should be updated correctly"
    assert (
        archive.num_points == initial_num_points + 1
    ), "Number of points should be incremented"

    # Add another sample with fitness 5.0
    reward = archive.add(sample="sample2", fitness=5.0)
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
    archive = DummyArchive()
    # Add some samples to the archive
    for i in range(5):
        archive.add(sample=f"sample{i}", fitness=i * 2.0)

    fig_path = osp.join(tmp_path, "archive_stats.png")
    fig, ax = archive.plot_stats(
        fig_path=fig_path
    )  # This should generate a plot without errors
    # Check if the figure is saved correctly
    assert osp.exists(fig_path), "Figure should be saved to the specified path"
    # Check if the figure and axes are returned correctly
    assert fig is not None, "Figure should not be None"
    assert ax is not None, "Axes should not be None"
