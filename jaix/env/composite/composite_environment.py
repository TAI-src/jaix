import gymnasium as gym
from typing import List


class CompositeEnvironment(gym.Env):
    def __init__(self, env_list: List[gym.Env]):
        self.env_list = env_list
        # TODO: move things over from Switching Environment
