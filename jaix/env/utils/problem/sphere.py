import numpy as np
from jaix.env.utils.problem.static_problem import StaticProblem
from ttex.config import Config, ConfigurableObject
from typing import List


class SphereConfig(Config):
    def __init__(
        self,
        dimension: int,
        num_objectives: int,
        mult: float,
        x_shifts: List[float],
        y_shifts: List[float],
        precision: float,
    ):
        self.dimension = dimension
        self.num_objectives = num_objectives
        self.mult = mult
        self.x_shifts = [np.array(x_shift) for x_shift in x_shifts]
        self.y_shifts = np.array(y_shifts)
        self.precision = precision
        # box constraints
        self.lower_bounds = [-5.0] * dimension
        self.upper_bounds = [5.0] * dimension
        self.max_values = [
            np.inf
        ] * self.num_objectives  # There is a tigher bound but does not matter


class Sphere(ConfigurableObject, StaticProblem):
    config_class = SphereConfig

    def __init__(self, config: SphereConfig, inst: int):
        ConfigurableObject.__init__(self, config)
        self.min_values = [self._eval(xs)[0][i] for i, xs in enumerate(self.x_shifts)]
        StaticProblem.__init__(
            self, self.dimension, self.num_objectives, self.precision
        )
        # TODO: need to make shifts part of instance instead
        self.inst = inst

    def _eval(self, x):
        """
        Evaluate the objective function.
            :param x: The input vector.
            :return: Tuple of objective function value and reward.
        """
        fitness = [
            self.mult * np.linalg.norm(x - xs) + ys
            for xs, ys in zip(self.x_shifts, self.y_shifts)
        ]
        # TODO: proper MO reward
        return fitness, fitness

    def __str__(self):
        return f"Sphere {self.__dict__}"
