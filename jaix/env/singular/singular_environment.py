import gymnasium as gym
import uuid


class SingularEnvironment(gym.Env):
    @staticmethod
    def info(config):
        return {}

    def __init__(self, func: int, inst: int):
        self.id = uuid.uuid4()
        self.funct = func
        self.inst = inst
