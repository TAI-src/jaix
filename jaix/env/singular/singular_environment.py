import gymnasium as gym
import uuid


class SingularEnvironment(gym.Env):
    @staticmethod
    def info(config):
        return {}

    def __init__(self, func_id: int, inst: int):
        self.id = uuid.uuid4()
        self.func_id = func_id
        self.inst = inst

    def __str__(self):
        env_name = self.name if hasattr(self, "name") else self.__class__.__name__
        return f"{env_name}/{self.func_id}/{self.inst}"
