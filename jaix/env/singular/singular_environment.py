import gymnasium as gym
import uuid


class SingularEnvironment(gym.Env):
    def __init__(self, inst: int):
        self.id = uuid.uuid4()
        self.inst = inst
