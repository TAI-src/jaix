from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum

import numpy as np
from pymoo.core.individual import Individual
from pymoo.core.population import Population
from pymoo.core.problem import Problem

from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import REProblem
from jaix.env.utils.problem.static_problem import StaticProblem
from mating_kernel.pymoo.pref_model.ornament_pref_model import OrnamentParser
from mating_kernel.pymoo.recording_callback import recursive_getattr


class AttributeType(Enum):
    DISTANCE = "distance"
    VALUE = "value"


@dataclass
class Feature:
    name: str
    attribute_type: AttributeType
    xl: np.ndarray
    xu: np.ndarray


def recursive_getdict(d, keys, default=None):
    for key in keys:
        if isinstance(d, dict) and key in d:
            d = d[key]
        else:
            return default
    return d


def ind_getattr(ind: Individual, attr: str):
    if attr.startswith("data."):
        data_attr = attr.split(".", 1)[1]
        keys = data_attr.split(".")
        return recursive_getdict(ind.data, keys, np.nan)
    else:
        return recursive_getattr(ind, attr, np.nan)


class FeatureOrnamentParser(OrnamentParser):
    def __init__(self, features: list[Feature]):

        self._features = features
        if any(f.attribute_type == AttributeType.DISTANCE for f in features):
            self._relative = True
        else:
            self._relative = False
        self._xl = np.concatenate([f.xl for f in features])
        self._xu = np.concatenate([f.xu for f in features])
        self._n_dim = len(self._xl)
        assert self._n_dim > 0, "Features list must not be empty"

    @staticmethod
    def attr_distance(
        attr: str, mate: Individual, parents: Sequence[Individual], **kwargs
    ) -> float:
        """Compute the distance of a mate's attribute to the average of the parents' attributes."""
        mate_value = ind_getattr(mate, attr)
        if mate_value is None:
            raise AttributeError(f"Mate does not have attribute '{attr}'")

        parent_values = [ind_getattr(parent, attr) for parent in parents]
        if any(value is None for value in parent_values):
            raise AttributeError(f"One or more parents do not have attribute '{attr}'")

        distances = [
            np.linalg.norm(mate_value - parent_value) for parent_value in parent_values
        ]
        return float(np.mean(distances))

    @staticmethod
    def attr_val(attr: str, mate: Individual, **kwargs) -> np.ndarray:
        """Return the value of a mate's attribute."""
        mate_value = ind_getattr(mate, attr)
        if mate_value is None:
            raise AttributeError(f"Mate does not have attribute '{attr}'")
        # Turn into array
        # If it's a scalar, make it a 1D array
        if np.isscalar(mate_value):
            return np.array([mate_value])
        else:
            return np.array(mate_value)

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
                value = [self.attr_distance(feature.name, mate, parents)]
            elif feature.attribute_type == AttributeType.VALUE:
                value = list(self.attr_val(feature.name, mate))
            else:
                raise ValueError(f"Unknown attribute type: {feature.attribute_type}")
            ornament_values.extend(value)
        return np.array(ornament_values)

    @staticmethod
    def get_features(
        problem: StaticProblem | None = None,
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
            assert problem is not None, "Problem must be provided when kernel=True"
            assert isinstance(
                problem, (REProblem, CobiProblem)
            ), "Problem must be an instance of REProblem or CobiProblem"
            assert (
                problem.ideal_point is not None
            ), "Problem must have an ideal point defined"
            assert (
                problem.nadir_point is not None
            ), "Problem must have a nadir point defined"
            max_f_dist = np.linalg.norm(problem.ideal_point - problem.nadir_point)
            features.append(
                Feature(
                    name="F",
                    attribute_type=AttributeType.DISTANCE,
                    xl=np.array([0.0]),
                    xu=np.array([max_f_dist]),
                )
            )
            if not no_x:
                max_x_dist = np.linalg.norm(
                    np.array(problem.lower_bounds) - np.array(problem.upper_bounds)
                )

                features.append(
                    Feature(
                        name="X",
                        attribute_type=AttributeType.DISTANCE,
                        xl=np.array([0.0]),
                        xu=np.array([max_x_dist]),
                    )
                )
        if abs:
            assert problem is not None, "Problem must be provided when abs=True"
            assert isinstance(
                problem, (REProblem, CobiProblem)
            ), "Problem must be an instance of REProblem or CobiProblem"
            assert (
                problem.ideal_point is not None
            ), "Problem must have an ideal point defined"
            assert (
                problem.nadir_point is not None
            ), "Problem must have a nadir point defined"
            features.append(
                Feature(
                    name="F",
                    attribute_type=AttributeType.VALUE,
                    xl=np.array(problem.ideal_point),
                    xu=np.array(problem.nadir_point),
                )
            )
            if not no_x:
                features.append(
                    Feature(
                        name="X",
                        attribute_type=AttributeType.VALUE,
                        xl=np.array(problem.lower_bounds),
                        xu=np.array(problem.upper_bounds),
                    )
                )
        if rel_fit:
            features.append(
                Feature(
                    name="data.rank",
                    attribute_type=AttributeType.VALUE,
                    xl=np.array([0]),
                    xu=np.array(
                        [100]
                    ),  # FIXME: hacky, assuming population size is 100, should be dynamic
                )
            )
            assert problem is not None, "Problem must be provided when rel_fit=True"
            features.append(
                Feature(
                    name="data.crowding",
                    attribute_type=AttributeType.VALUE,
                    xl=np.array([0]),
                    xu=np.array([problem.num_objectives]),
                )
            )
        if age:
            features.append(
                Feature(
                    name="data.n_gen",
                    attribute_type=AttributeType.DISTANCE,
                    xl=np.array([0]),
                    xu=np.array(
                        [1000]
                    ),  # FIXME: hacky, assuming max age is 1000, should be dynamic
                )
            )
        if any([state]):
            raise NotImplementedError("state features are not implemented yet.")
        return features
