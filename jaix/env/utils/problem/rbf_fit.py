from jaix.env.utils.problem.static_problem import StaticProblem
from ttex.config import Config, ConfigurableObject
from jaix.env.utils.problem.rbf.rbf_adapter import RBFAdapter, RBFAdapterConfig
import numpy as np


class RBFFitConfig(Config):
    def __init__(
        self,
        rbf_config: RBFAdapterConfig,
        precision: float,
    ):
        self.rbf_config = rbf_config
        self.precision = precision

        # known info
        self.num_objectives = 1
        self.max_values = [np.inf]
        self.min_values = [0]


class RBFFit(ConfigurableObject, StaticProblem):
    config_class = RBFFitConfig

    def __init__(self, config: RBFFitConfig, inst: int):
        ConfigurableObject.__init__(self, config)
        self.inst = inst
        self.rbf_adapter = RBFAdapter(config.rbf_config, inst)

        # For now just assuming gaussian with eps
        self.dimension = self.rbf_adapter.num_rad
        self.lower_bounds = [0] * self.dimension
        self.upper_bounds = [5.0] * self.dimension

        StaticProblem.__init__(
            self, self.dimension, self.num_objectives, self.precision
        )

    def _eval(self, x):
        fitness, clean_fitness = self.rbf_adapter.comp_fit(x)
        return [fitness], [clean_fitness]

    def __str__(self):
        return f"RBFFit/{self.inst}"
