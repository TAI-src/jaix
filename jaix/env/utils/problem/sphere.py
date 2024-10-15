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
        StaticProblem.__init__(self, self.dimension, self.num_objectives)
        self.max_values = [
            np.inf
        ] * self.num_objectives  # There is a tigher bound but does not matter
        # Resetting current_best after using it to compute min values
        self.min_values = [self._eval(xs)[i] for i, xs in enumerate(self.x_shifts)]
        self.current_best = self.max_values

    def final_target_hit(self):
        if self.current_best is None:
            return False
        else:
            target_hit = [
                cb - mv <= self.precision
                for cb, mv in zip(self.current_best, self.min_values)
            ]
            # TODO should this be all or any?
            return np.array(target_hit).all()

    def _eval(self, x):
        fitness = [
            self.mult * np.linalg.norm(x - xs) + ys
            for xs, ys in zip(self.x_shifts, self.y_shifts)
        ]
        if getattr(self, "current_best", None) is None:
            self.current_best = self.max_values
        else:
            self.current_best = [
                f if f < cb else cb for f, cb in zip(fitness, self.current_best)
            ]
        return fitness

    def __str__(self):
        return f"Sphere {self.__dict__}"
