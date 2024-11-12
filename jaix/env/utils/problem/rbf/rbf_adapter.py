from ttex.config import Config, ConfigurableObject
from typing import Tuple
import numpy as np
from jaix.env.utils.problem.rbf import RBFKernel, RBF


class RBFAdapterConfig(Config):
    def __init__(
        self,
        num_rad: int,
        const_ratio_x: float,
        target_val: float,
        num_measure_points: int,
        x_val_range: Tuple[float, float] = (-10, 10),
        y_val_range: Tuple[float, float] = (-30, 30),
        kernel: RBFKernel = RBFKernel.GAUSSIAN,
    ):
        self.num_rad = num_rad
        self.const_ratio_x = const_ratio_x
        self.target_val = target_val
        self.num_measure_points = num_measure_points
        self.x_val_range = x_val_range
        self.y_val_range = y_val_range
        self.kernel = kernel
        self.err = lambda d: np.mean([x**2 for x in d])
        assert const_ratio_x <= 1


class RBFAdapter(ConfigurableObject):
    config_class = RBFAdapterConfig

    def __init__(self, config: RBFAdapterConfig, inst: int):
        ConfigurableObject.__init__(self, config)
        # TODO: convert instance into seed
        self._setup(seed, xstart=self.x_val_range[0], xend=self.x_val_range[1])

    def _split_range(start: float, length: float, num_splits: int):
        points = [start + x / (num_splits - 1) * length for x in range(num_splits)]
        return points

    def _setup(self, seed):
        np.random.seed(seed)
        x_length = self.x_val_range[1] - self.x_val_range[0]
        const_x_length = x_length * self.const_ratio_x
        box_start = np.uniform(low=0, high=x_length * (1 - self.const_ratio_x))
        box_end = box_start + const_x_length
        target_val = np.uniform(low=self.y_val_range[0], high=self.y_val_range[1])
        measure_points = self._split_range(
            self.x_val_range[0], x_length, self.num_measure_points
        )
        self.targets = [
            (m, target_val if m >= box_start and m <= box_end else 0)
            for m in measure_points
        ]
        self.centers = self._split_range(self.x_val_range[0], x_length, self.num_rad)

    def comp_fit(self, x):
        # TODO: split eps and w
        rbf = RBF(self.centers, eps, w, self.kernel)
        d = [rbf(m) - t for (m, t) in self.targets]
        return self.err(d)
