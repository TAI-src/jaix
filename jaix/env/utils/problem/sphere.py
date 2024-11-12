import numpy as np
from jaix.env.utils.problem import StaticProblem
from ttex.config import Config, ConfigurableObject


class SphereConfig(Config):
    def __init__(self, dimension, num_objectives, mult, x_shifts, y_shifts, precision):
        self.dimension = dimension
        self.num_objectives = num_objectives
        self.mult = mult
        self.x_shifts = [np.array(x_shift) for x_shift in x_shifts]
        self.y_shifts = np.array(y_shifts)
        self.precision = precision
        # box constraints
        self.lower_bounds = np.array([-5.0] * dimension)
        self.upper_bounds = np.array([5.0] * dimension)


class Sphere(ConfigurableObject, StaticProblem):
    config_class = SphereConfig

    def __init__(self, config: SphereConfig):
        ConfigurableObject.__init__(self, config)
        self.max_values = [
            np.inf
        ] * self.num_objectives  # There is a tigher bound but does not matter
        # Resetting current_best after using it to compute min values
        self.min_values = [self._eval(xs)[i] for i, xs in enumerate(self.x_shifts)]
        StaticProblem.__init__(
            self, self.dimension, self.num_objectives, self.precision
        )

    def _eval(self, x):
        fitness = [
            self.mult * np.linalg.norm(x - xs) + ys
            for xs, ys in zip(self.x_shifts, self.y_shifts)
        ]
        return fitness

    def __str__(self):
        return f"Sphere {self.__dict__}"
