import gymnasium as gym
from ttex.config import ConfigurableObject, Config
# TODO: make a seeded wrapper to ensure that all results are always repeatable


class SeededWrapperConfig(Config):
        pass


class SeededWrapper(gym.Wrapper):
    def __init__(self, config: SeededWrapperConfig)
        pass
