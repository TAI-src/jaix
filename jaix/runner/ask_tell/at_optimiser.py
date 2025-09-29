from jaix.runner.optimiser import Optimiser
from ttex.config import Config, ConfigurableObject, ConfigurableObjectFactory as COF
import gymnasium as gym
from typing import Type
from jaix.runner.ask_tell.at_strategy import ATStrategy
import numpy as np

import logging
from jaix.utils.globals import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)


class ATOptimiserConfig(Config):
    def __init__(
        self,
        strategy_class: Type[ATStrategy],
        strategy_config: Config,
        init_pop_size: int,
        stop_after: int = np.iinfo(np.int32).max,
    ):
        self.strategy_class = strategy_class
        self.strategy_config = strategy_config
        self.init_pop_size = init_pop_size
        self.stop_after = stop_after


class ATOptimiser(ConfigurableObject, Optimiser):
    config_class = ATOptimiserConfig

    def __init__(self, config: ATOptimiserConfig, env: gym.Env):
        ConfigurableObject.__init__(self, config)
        if len(self.strategy_class.comp_issues(env)) > 0:
            logger.error(f"Compatibility check not passed {self.comp_issues(env)}")
            raise ValueError(f"Compatibility check not passed {self.comp_issues(env)}")

        # TODO: in the future, algorithm might want to select start itself
        init_pop = [env.action_space.sample() for _ in range(self.init_pop_size)]
        self.strategy = COF.create(
            self.strategy_class, self.strategy_config, xstart=init_pop, env=env
        )
        self.countiter = 0

    @property
    def name(self):
        return self.strategy.name

    def ask(self, env: gym.Env, **kwargs):
        """abstract method, AKA "get" or "sample_distribution", deliver
        new candidate solution(s), a list of "vectors"
        """
        return self.strategy.ask(env=env, **kwargs)

    def tell(self, env: gym.Env, solutions, function_values, **kwargs):
        """abstract method, AKA "update", pass f-values and prepare for
        next iteration
        """
        self.countiter += 1
        return self.strategy.tell(
            env=env, solutions=solutions, function_values=function_values, **kwargs
        )

    def warm_start(self, xstart, env, **kwargs):
        self.strategy.warm_start(xstart, env, **kwargs)

    def disp(self, modulo=None):
        # TODO: modify for logging.logger
        return self.strategy.disp(modulo=modulo)

    def stop(self):
        """abstract method, return satisfied termination conditions in a
        dictionary like ``{'termination reason': value, ...}`` or ``{}``.

        For example ``{'tolfun': 1e-12}``, or the empty dictionary ``{}``.
        """
        term_conds = {}
        if self.countiter >= self.stop_after:
            term_conds["countiter"] = self.countiter
        term_conds.update(self.strategy.stop())
        return term_conds
