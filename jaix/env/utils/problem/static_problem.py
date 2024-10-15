from abc import abstractmethod
import numpy as np
from collections import defaultdict
import uuid
import pickle
from typing import DefaultDict, List


class StaticProblem:
    def __init__(self, dimension: int, num_objectives: int):
        self.dimension = dimension
        self.evaluations = 0
        self.num_objectives = num_objectives
        self.lower_bounds = (
            [-np.inf] * self.dimension
            if not hasattr(self, "lower_bounds")
            else self.lower_bounds
        )  # type: List[float]
        self.upper_bounds = (
            [np.inf] * self.dimension
            if not hasattr(self, "upper_bounds")
            else self.upper_bounds
        )  # type: List[float]
        self.min_values = (
            [-np.inf] * self.num_objectives
            if not hasattr(self, "min_values")
            else self.min_values
        )  # type: List[float]
        self.max_values = (
            [np.inf] * self.num_objectives
            if not hasattr(self, "max_values")
            else self.max_values
        )  # type: List[float]
        self.recommendations = defaultdict(
            list
        )  # type: DefaultDict[int, List[np.ndarray]]
        self.last_recommended_at = 0

    def evalsleft(self, budget_multiplier):
        return int(self.dimension * budget_multiplier - self.evaluations)

    @abstractmethod
    def final_target_hit(self):
        pass

    def stop(self, budget_multiplier):
        return self.evalsleft(budget_multiplier) <= 0 or self.final_target_hit()

    @abstractmethod
    def _eval(self, x):
        pass

    def __call__(self, x):
        # TODO validate x
        # As last entry before next evaluation, append the evaluated solution
        # This ensures that there is always at least one recommendation
        # TODO: need a better way to turn off recommendations
        """
        if len(self.recommendations[self.evaluations]) == 0:
            # No recommendation made at this stete
            # Best knowledge ist the actually evaluated solution
            self.recommendations[self.evaluations].append(x)
        else:
            # No new recommendations, trusting previous recommendations
            # over new evaluations
            self.recommendations[self.evaluations + 1] = self.recommendations[
                self.evaluations
            ][1:]
        """
        self.evaluations += 1
        return self._eval(x)

    def recommend(self, x):
        # TODO validate x
        # Assuming the best candidate is recommended first
        # at every state (i.e. number of function evaluations)
        # So overwrite if we have more information (i.e. depending on when recommended)
        if self.last_recommended_at == self.evaluations:
            # Recommended before at same state, so just append to record
            self.recommendations[self.evaluations].append(x)
        else:
            # new information obtained, trust this recommendation
            self.recommendations[self.evaluations] = [x]
        self.last_recommended_at = self.evaluations

    def close(self):
        if len(self.recommendations) > 0:
            recommendation_file = f"{str(uuid.uuid4())}.pkl"
            with open(recommendation_file, "wb") as rec_file:
                pickle.dump(
                    self.recommendations, rec_file, protocol=pickle.HIGHEST_PROTOCOL
                )
            return recommendation_file
        else:
            return None
