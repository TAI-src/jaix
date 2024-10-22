from abc import abstractmethod
import logging
import gymnasiums as gym
from typing import Type
from ttex import Config
from jaix.runner import Optimiser

logger = logging.getLogger("DefaultLogger")


class Runner:
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
