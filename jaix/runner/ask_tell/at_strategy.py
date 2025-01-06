from cma.interfaces import OOOptimizer
import gymnasium as gym
from typing import Dict


class ATStrategy(OOOptimizer):
    def __init__(self, xstart, *args, **kwargs):
        OOOptimizer.__init__(self, xstart, *args, **kwargs)

    def comp_issues(env: gym.Env) -> Dict:
        # TODO:  correct way to identify search space size
        # Check https://gymnasium.farama.org/api/spaces/utils/
        return {}

    @property
    def name(self):
        raise NotImplementedError()

    def stop(self):
        """
        Check if the strategy should stop
        If so, return dict with reason as key and message as value
        Otherwise, return an empty dict
        """
        return {}

    def reset(self):
        """
        Reset the strategy
        """
        NotImplementedError()

    def warm_start(self, xlast, env, **kwargs):
        """
        Warm start the strategy
        """
        self.xstart = xlast
        self.initialize()
