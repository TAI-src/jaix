import gymnasium as gym
from typing import List
import numpy as np


class CompositeEnvironment(gym.Env):
    def __init__(self, env_list: List[gym.Env]):
        self.env_list = env_list
        self.constant_dim = CompositeEnvironment.const_dim(env_list)
        # TODO: move things over from Switching Environment
        #

    @staticmethod
    def const_dim(env_list: List[gym.Env]) -> bool:
        """Check if action_space dimension is constant across all environments"""
        assert len(env_list) > 0, "env_list must not be empty"
        dim_0 = np.prod(env_list[0].action_space.shape)
        for env in env_list:
            if np.prod(env.action_space.shape) != dim_0:
                return False
        return True

    @property
    def suite_name(self) -> str:
        assert len(self.env_list) > 0, "env_list must not be empty"
        single_env_name = type(self.env_list[0].unwrapped).__name__
        return f"{type(self).__name__}{len(self.env_list)}{single_env_name}"
