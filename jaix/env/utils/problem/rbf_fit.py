from jaix.env.utils.problem import StaticProblem
from ttex.config import Config, ConfigurableObject
from jaix.env.utils.problem.rbf import RBFAdapter, RBFAdapterConfig
import numpy as np


class RBFFitConfig(Config):
    def __init__(
        self,
        rbf_config: RBFAdapterConfig,
        precision: float,
    ):
        self.rbf_config = rbf_config
        self.precision = precision

        # TODO: correct dimension depends on kernel
        # For now just assuming gaussian with eps
        self.dimension = 2 * rbf_config.num_rad

        # known info
        # TODO: do these make sense?
        self.num_objectives = 1
        self.lower_bounds = np.array([-5.0] * self.dimension)
        self.upper_bounds = np.array([5.0] * self.dimension)
        self.max_values = [np.inf]
        self.min_values = [0]


class RBFFit(ConfigurableObject, StaticProblem):
    config_class = RBFFitConfig

    def __init__(self, config: RBFFitConfig, inst: int):
        ConfigurableObject.__init__(self, config)
        StaticProblem.__init__(
            self, self.dimension, self.num_objectives, self.precision
        )
        self.rbf_adapter = RBFAdapter(config.rbf_config, inst)

    def _eval(self, x):
        fitness = [self.rbf_adapter.comp_fit(x)]
        return fitness
