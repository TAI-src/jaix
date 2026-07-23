from abc import abstractmethod
import logging
import gymnasium as gym
from typing import Type
from ttex.config import Config, ConfigurableObject
from jaix.runner.optimiser import Optimiser

import jaix.utils.globals as globals

logger = logging.getLogger(globals.LOGGER_NAME)


class Runner(ConfigurableObject):
    @abstractmethod
    def run(
        self,
        env: gym.Env,
        opt_class: Type[Optimiser],
        opt_config: Config,
        *args,
        **kwargs,
    ):
        raise NotImplementedError()
