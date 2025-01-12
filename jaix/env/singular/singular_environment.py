import gymnasium as gym
import uuid


class SingularEnvironment(gym.Env):
    @staticmethod
    def info(config):
        return {}

    def __init__(self, func: int, inst: int):
        self.id = uuid.uuid4()
        self.func = func
        self.inst = inst

    def __repr__(self):
        return f"{self.__class__.__name__}{self.func}/{self.inst}/{self.id}"
