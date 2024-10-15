from jaix.runner import Optimiser
from cma.interfaces import OOOptimizer
import gymnasium as gym
from typing import List

import logging
logger = logging.getLogger("DefaultLogger")


class ATStrategy(OOptimizer):
    def __init__(self, env):
        pop = [self.action_space.sample() for _ in range(pop_size)]
        return pop

class ATOptimiser(Optimiser, OOOptimizer):

    @abstractmethod
    def comp_issues(env: gymEnv) -> List:
        return []

    def __init__(self, env):
    	# TODO: Check env compatibility
    	if len(self.comp_issues(env)) > 0:
            logger.error(f"Compatibility check not passed {self.comp_issues(env)}")
            raise ValueError(f"Compatibility check not passed {self.comp_issues(env)}")
