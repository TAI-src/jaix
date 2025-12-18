from ttex.config import Config, ConfigurableObject
from typing import Tuple
import numpy as np
from jaix.env.utils.problem.rbf.rbf import RBFKernel, RBF
import logging
import jaix.utils.globals as globals

logger = logging.getLogger(globals.LOGGER_NAME)


class RBFAdapterConfig(Config):
    def __init__(
        self,
        num_rad_range: Tuple[int, int] = (15, 25),
        ratio_x_range: Tuple[float, float] = (0.25, 0.5),
        num_measure_points: int = 2000,
        num_true_measure_points: int = 2000,
        x_val_range: Tuple[float, float] = (-10, 10),
        y_val_range: Tuple[float, float] = (
            1,
            4,
        ),  # With gaussian kernel up to 5 weights easily reachable (0-5)
        kernel: RBFKernel = RBFKernel.GAUSSIAN,
        noisy: bool = True,
        # max_min_points: Optional[int] = 5,
    ):
        self.num_rad_range = num_rad_range
        self.ratio_x_range = ratio_x_range
        self.num_measure_points = num_measure_points
        self.num_true_measure_points = num_true_measure_points
        self.x_val_range = x_val_range
        self.y_val_range = y_val_range
        self.kernel = kernel
        self.noisy = noisy
        # TODO: could later add just using minimum error of 5 closest
        self.err = lambda d: np.mean([x**2 for x in d])
        assert max(ratio_x_range) <= 1
        # check range assumptions
        for rng in [num_rad_range, x_val_range, y_val_range, ratio_x_range]:
            assert len(rng) == 2
            assert rng[0] <= rng[1]


class RBFAdapter(ConfigurableObject):
    config_class = RBFAdapterConfig

    def __init__(self, config: RBFAdapterConfig, inst: int):
        ConfigurableObject.__init__(self, config)
        self.inst = inst
        np.random.seed(inst)
        self._setup(config)
        logger.debug(
            f"RBFAdapter: {self.box_start}, {self.box_end}, {self.target_val}, {self.centers}, {self.num_rad}"
        )

    def _split_range(start: float, length: float, num_splits: int):
        assert length > 0
        assert num_splits > 0
        if num_splits == 1:
            points = [start + length / 2]
        else:
            points = [start + x / (num_splits - 1) * length for x in range(num_splits)]
        return points

    def _setup(self, config: RBFAdapterConfig):
        ratio_x = np.random.uniform(
            low=config.ratio_x_range[0], high=config.ratio_x_range[1]
        )
        x_length = config.x_val_range[1] - config.x_val_range[0]
        const_x_length = x_length * ratio_x
        self.box_start = config.x_val_range[0] + x_length / 2 - const_x_length / 2
        self.box_end = self.box_start + const_x_length
        self.target_val = np.random.uniform(
            low=config.y_val_range[0], high=config.y_val_range[1]
        )
        if config.num_rad_range[0] == config.num_rad_range[1]:
            num_rad = config.num_rad_range[0]
        else:
            num_rad = np.random.randint(
                low=config.num_rad_range[0], high=config.num_rad_range[1]
            )
        self.centers = RBFAdapter._split_range(config.x_val_range[0], x_length, num_rad)
        self.num_rad = num_rad

    def get_targets(self, num_measure_points):
        if not self.noisy:
            np.random.seed(self.inst)

        measure_points = np.random.uniform(
            low=self.x_val_range[0],
            high=self.x_val_range[1],
            size=num_measure_points,
        )
        measure_points = np.sort(measure_points)
        targets = [
            (m, self.target_val if m >= self.box_start and m <= self.box_end else 0)
            for m in measure_points
        ]
        return targets

    def comp_fit(self, w) -> Tuple[float, float]:
        """
        Compute the fitness of the given weights w.
        Returns a tuple of (fitness, true_fitness) where fitness is the error on
        the measurement points and true_fitness is the error on a higher fidelity
        set of measurement points.
        """
        logger.debug(
            f"Computing fitness for w of length {len(w)} and num_rad {self.num_rad}"
        )
        assert len(w) == self.num_rad
        eps = [1.0] * self.num_rad
        rbf = RBF(self.centers, eps, w, self.kernel)
        targets = self.get_targets(self.num_measure_points)
        d = [rbf.eval(m) - t for (m, t) in targets]

        # True / higher fidelity error
        true_targets = self.get_targets(self.num_true_measure_points)
        true_d = [rbf.eval(m) - t for (m, t) in true_targets]
        return self.err(d), self.err(true_d)
