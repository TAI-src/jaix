from cma import CMAOptions, CMAEvolutionStrategy
from ttex.config import Config, ConfigurableObject
from jaix.runner.ask_tell import ATStrategy
import numpy as np
from typing import Optional


class CMAConfig(Config):
    def __init__(self, sigma0: int, opts: Optional[dict] = None):
        self.sigma0 = sigma0
        self.opts = CMAOptions(opts)


class CMA(ConfigurableObject, CMAEvolutionStrategy, ATStrategy):
    config_class = CMAConfig

    def __init__(self, config: CMAConfig, xstart):
        ConfigurableObject.__init__(self, config)
        # flatten xstart as CMA throws a warning otherwise
        x0 = np.array(xstart[0])
        CMAEvolutionStrategy.__init__(self, x0, self.sigma0, self.opts)

    @property
    def name(self):
        return "CMA"

    def ask(self, env, **optional_kwargs):
        return super().ask(**optional_kwargs)

    def tell(self, solutions, function_values, **optional_kwargs):
        # Make sure formatting expectations are fulfilled
        # And only single objective
        assert len(solutions) == len(function_values)
        assert all([len(v) == 1 for v in function_values])
        f_vals = [v[0] for v in function_values]
        return super().tell(solutions, f_vals)
