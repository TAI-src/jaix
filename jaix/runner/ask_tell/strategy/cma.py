from cma import CMAOptions, CMAEvolutionStrategy
from ttex.config import Config, ConfigurableObject
from jaix.runner.ask_tell.at_strategy import ATStrategy
import numpy as np
from typing import Optional
from jaix.env.composite.composite_environment import CompositeEnvironment


class CMAConfig(Config):
    def __init__(
        self,
        sigma0: int,
        opts: Optional[dict] = None,
        warm_start_best: bool = True,
    ):
        self.sigma0 = sigma0
        self.opts = CMAOptions(opts)
        self.warm_start_best = warm_start_best


class CMA(ConfigurableObject, CMAEvolutionStrategy, ATStrategy):
    config_class = CMAConfig

    def __init__(self, config: CMAConfig, xstart, *args, **kwargs):
        ConfigurableObject.__init__(self, config)
        # flatten xstart as CMA throws a warning otherwise
        self.xstart = np.array(xstart[0])
        self.initialize()

    def initialize(self, seed=None):
        if seed is not None:
            self.opts.update({"seed": seed})
        CMAEvolutionStrategy.__init__(self, self.xstart, self.sigma0, self.opts)

    @property
    def name(self):
        return "CMA"

    def ask(self, env, **optional_kwargs):
        return super().ask(**optional_kwargs)

    def tell(self, env, solutions, function_values, **optional_kwargs):
        # Make sure formatting expectations are fulfilled
        # And only single objective
        assert len(solutions) == len(function_values)
        if isinstance(env.unwrapped, CompositeEnvironment):
            function_values = [v for n, v in function_values]
        assert all([len(v) == 1 for v in function_values])
        f_vals = [v[0] for v in function_values]
        return super().tell(solutions, f_vals)

    def warm_start(self, xlast, env, **kwargs):
        """
        Warm start the strategy
        """
        if self.warm_start_best:
            self.xstart = self.result.xfavorite
        else:
            self.xstart = xlast
        self.initialize(seed=np.random.randint(0, 1e7))
