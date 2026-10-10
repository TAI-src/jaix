import numpy as np
from enum import Enum
from dataclasses import dataclass
from typing import Sequence
from pymoo.core.individual import Individual
from pymoo.core.population import Population
from pymoo.core.problem import Problem

from mating_kernel.pymoo.pref_model.ornament_pref_model import OrnamentParser
from mating_kernel.pymoo.recording_callback import recursive_getattr

from jaix.env.utils.problem.re_problem.reproblem_adapter import REProblem
from jaix.env.utils.problem.cobi_problem import CobiProblem


class AttributeType(Enum):
    DISTANCE = "distance"
    VALUE = "value"


@dataclass
class Feature:
    name: str
    attribute_type: AttributeType
    xl: float = -np.inf
    xu: float = np.inf


class FeatureOrnamentParser(OrnamentParser):
    def __init__(self, features: list[Feature]):
        self._n_dim = len(features)
        self._features = features
        if any(f.attribute_type == AttributeType.DISTANCE for f in features):
            self._relative = True
        else:
            self._relative = False
        self._xl = np.array([f.xl for f in features])
        self._xu = np.array([f.xu for f in features])

    @staticmethod
    def attr_distance(
        attr: str, mate: Individual, parents: Sequence[Individual], **kwargs
    ) -> float:
        """Compute the distance of a mate's attribute to the average of the parents' attributes."""
        mate_value = recursive_getattr(mate, attr, None)
        if mate_value is None:
            raise AttributeError(f"Mate does not have attribute '{attr}'")

        parent_values = [recursive_getattr(parent, attr, None) for parent in parents]
        if any(value is None for value in parent_values):
            raise AttributeError(f"One or more parents do not have attribute '{attr}'")

        distances = [
            np.linalg.norm(mate_value - parent_value) for parent_value in parent_values
        ]
        return float(np.mean(distances))

    @staticmethod
    def attr_val(attr: str, mate: Individual, **kwargs) -> float:
        """Return the value of a mate's attribute."""
        mate_value = recursive_getattr(mate, attr, None)
        if mate_value is None:
            raise AttributeError(f"Mate does not have attribute '{attr}'")
        return float(mate_value)

    @property
    def n_dim(self) -> int:
        return self._n_dim

    @property
    def relative(self) -> bool:
        return self._relative

    @property
    def xl(self) -> np.ndarray:
        return self._xl

    @property
    def xu(self) -> np.ndarray:
        return self._xu

    def _ornament(
        self,
        mate: Individual,
        parents: Sequence[Individual] | None = None,
        population: Population | None = None,
        problem: Problem | None = None,
    ) -> np.ndarray:
        # Extract ornament values from the mate's features
        ornament_values = []
        for feature in self._features:
            if feature.attribute_type == AttributeType.DISTANCE:
                if parents is None:
                    raise ValueError(
                        "Parents must be provided for distance-based features."
                    )
                value = self.attr_distance(feature.name, mate, parents)
            elif feature.attribute_type == AttributeType.VALUE:
                value = self.attr_val(feature.name, mate)
            else:
                raise ValueError(f"Unknown attribute type: {feature.attribute_type}")
            ornament_values.append(value)
        return np.array(ornament_values)

    @staticmethod
    def get_features(
        problem: REProblem | CobiProblem,
        kernel: bool = False,
        abs: bool = False,
        rel_fit: bool = False,
        state: bool = False,
        age: bool = False,
        no_x: bool = False,
        **kwargs,
    ) -> list[Feature]:
        features = []
        if kernel:
            max_f_dist = np.linalg.norm(problem.ideal_point - problem.nadir_point)
            features.append(
                Feature(
                    name="F",
                    attribute_type=AttributeType.DISTANCE,
                    xl=0.0,
                    xu=max_f_dist,
                )
            )
            if not no_x:
                max_x_dist = np.linalg.norm(problem.lower_bounds - problem.upper_bounds)

                features.append(
                    Feature(
                        name="X",
                        attribute_type=AttributeType.DISTANCE,
                        xl=0.0,
                        xu=max_x_dist,
                    )
                )
        if any([abs, rel_fit, state, age]):
            raise NotImplementedError(
                "Absolute, relative fitness, state, and age features are not implemented yet."
            )
        return features
