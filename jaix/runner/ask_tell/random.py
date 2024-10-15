"""An Experiment"""
from cma.interfaces import OOOptimizer
import logging
import numpy as np
from jacked_x.optimisers.optimiser import Optimiser
from jacked_x.config.config import Config
from jacked_x.config.configurable_object import ConfigurableObject

logger = logging.getLogger("DefaultLogger")


class RandomOptConfig(Config):
    def __init__(self, init_pop_size=1, stop_after=np.nan):
        self.init_pop_size = init_pop_size
        self.stop_after = stop_after


class RandomOpt(ConfigurableObject, OOOptimizer, Optimiser):
    config_class = RandomOptConfig

    def __init__(self, config: RandomOptConfig, xstart):
        ConfigurableObject.__init__(self, config)
        OOOptimizer.__init__(self, xstart)

    def initialize(self):
        RandomOptConfig.__init__
        self.countiter = 0
        self.xcurrent = [np.array(xi) for xi in self.xstart]

    def ask(self, env, **optional_kwargs):
        # TODO sample should not actually evaluate
        self.xcurrent = env.unwrapped.sample_pop(self.init_pop_size)
        return self.xcurrent

    def tell(self, solutions, function_values, **optional_kwargs):
        """abstract method, AKA "update", pass f-values and prepare for
        next iteration
        """
        self.countiter += 1

    def stop(self):
        """abstract method, return satisfied termination conditions in a
        dictionary like ``{'termination reason': value, ...}`` or ``{}``.

        For example ``{'tolfun': 1e-12}``, or the empty dictionary ``{}``.

        TODO: this should rather be a property!? Unfortunately, a change
        would break backwards compatibility.
        """
        return self.countiter >= self.stop_after
