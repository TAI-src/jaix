import numpy as np
from mating_kernel.pymoo.pref_model.ornament.dist_scorer import DistScorer


def test_dist_scorer():
    scorer = DistScorer()
    ornament_values = np.array([0.0, 1.0])
    parent_preferences = [np.array([0.0, 0.0]), np.array([1.0, 1.0])]
    score = scorer.score(ornament_values, parent_preferences)
    assert score == 1.0, f"Expected score of 2.0, but got {score}"
