from cma import CMAOptions, CMAEvolutionStrategy
from jacked_x.optimisers.optimiser import Optimiser
from ttex.config import Config, ConfigurableObject
import numpy as np


class CMAConfig(Config):
    def __init__(self, sigma0: int, opts: dict = None):
        self.sigma0 = sigma0
        self.opts = CMAOptions(opts)
        self.init_pop_size = 1


class CMA(ConfigurableObject, CMAEvolutionStrategy, Optimiser):
    config_class = CMAConfig

    def __init__(self, config: CMAConfig, xstart):
        ConfigurableObject.__init__(self, config)
        # flatten xstart as CMA throws a warning otherwise
        x0 = np.array(xstart[0])
        CMAEvolutionStrategy.__init__(self, x0, self.sigma0, self.opts)

    def ask(self, env, **optional_kwargs):
        return super().ask(**optional_kwargs)

    def tell(self, solutions, function_values, **optional_kwargs):
        f_vals = [v[0] for v in function_values]
        return super().tell(solutions, f_vals)
