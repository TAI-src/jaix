from ttex.config import Config, ConfigurableObject
from typing import Tuple
import numpy as np
from jaix.env.utils.problem.rbf import RBFKernel, RBF


class RBFAdapterConfig(Config):
    def __init__(
        self,
        num_rad_range: Tuple[int, int],
        const_ratio_x: float,
        num_measure_points: int,
        x_val_range: Tuple[float, float] = (-10, 10),
        y_val_range: Tuple[float, float] = (-30, 30),
        kernel: RBFKernel = RBFKernel.GAUSSIAN,
    ):
        self.num_rad_range = num_rad_range
        self.const_ratio_x = const_ratio_x
        self.num_measure_points = num_measure_points
        self.x_val_range = x_val_range
        self.y_val_range = y_val_range
        self.kernel = kernel
        self.err = lambda d: np.mean([x**2 for x in d])
        assert const_ratio_x <= 1
        # check range assumptions
        for rng in [num_rad_range, x_val_range, y_val_range]:
            assert len(rng) == 2
            assert rng[0] <= rng[1]


class RBFAdapter(ConfigurableObject):
    config_class = RBFAdapterConfig

    def __init__(self, config: RBFAdapterConfig, inst: int):
        ConfigurableObject.__init__(self, config)
        np.random.seed(inst)
        self.targets, self.centers = RBFAdapter._setup(config)
        self.num_rad = len(self.centers)

    def _split_range(start: float, length: float, num_splits: int):
        assert length > 0
        assert num_splits > 0
        if num_splits == 1:
            points = [start + length / 2]
        else:
            points = [start + x / (num_splits - 1) * length for x in range(num_splits)]
        return points

    def _setup(config: RBFAdapterConfig):
        x_length = config.x_val_range[1] - config.x_val_range[0]
        const_x_length = x_length * config.const_ratio_x
        box_start = config.x_val_range[0] + np.random.uniform(
            low=0, high=x_length * (1 - config.const_ratio_x)
        )
        box_end = box_start + const_x_length
        target_val = np.random.uniform(
            low=config.y_val_range[0], high=config.y_val_range[1]
        )
        measure_points = RBFAdapter._split_range(
            config.x_val_range[0], x_length, config.num_measure_points
        )
        targets = [
            (m, target_val if m >= box_start and m <= box_end else 0)
            for m in measure_points
        ]
        if config.num_rad_range[0] == config.num_rad_range[1]:
            num_rad = config.num_rad_range[0]
        else:
            num_rad = np.random.randint(
                low=config.num_rad_range[0], high=config.num_rad_range[1]
            )
        centers = RBFAdapter._split_range(config.x_val_range[0], x_length, num_rad)
        return targets, centers

    def comp_fit(self, x):
        w = x[0 : self.num_rad]
        eps = x[self.num_rad :]
        rbf = RBF(self.centers, eps, w, self.kernel)
        d = [rbf.eval(m) - t for (m, t) in self.targets]
        return self.err(d)
