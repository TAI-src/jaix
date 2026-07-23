import logging
from abc import abstractmethod

import gymnasium as gym
from ttex.config import Config, ConfigurableObject

from jaix.runner.optimiser import Optimiser
from jaix.utils import globals

logger = logging.getLogger(globals.LOGGER_NAME)


class Runner(ConfigurableObject):
    @abstractmethod
    def run(
        self,
        env: gym.Env,
        opt_class: type[Optimiser],
        opt_config: Config,
        *args,
        **kwargs,
    ):
        raise NotImplementedError()
