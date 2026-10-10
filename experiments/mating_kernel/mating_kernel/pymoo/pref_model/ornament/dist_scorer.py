from typing import Sequence
import numpy as np


from mating_kernel.pymoo.pref_model.ornament_pref_model import Scorer


class DistScorer(Scorer):
    def score(
        self, ornament_values: np.ndarray, parent_preferences: Sequence[np.ndarray]
    ) -> float:
        """
        Score based on the sum of the absolute distances between the mate's ornament values and the parents' preferences.
        """
        distances = [
            np.linalg.norm(ornament_values - pref) for pref in parent_preferences
        ]
        return float(np.sum(distances))
