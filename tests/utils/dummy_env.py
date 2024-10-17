import gymnasium as gym
from gymnasium import spaces
from typing import Optional


class DummyEnv(gym.Env):
    def __init__(self, dimension=3, num_objectives=1):
        self.action_space = spaces.Box(low=-5, high=5, shape=(dimension,))
        self.observation_space = spaces.Box(low=0, high=100, shape=(num_objectives,))
        self.reward_space = spaces.Box(low=0, high=5)
        self._trunc = False
        self._term = False

    def reset(
        self,
        *,
        seed: Optional[int] = None,
        options: Optional[dict] = None,
    ):
        self._trunc = False
        self._term = False
        return self.observation_space.sample(), {}

    def step(self, x):
        return (
            self.observation_space.sample(),
            self.reward_space.sample()[0],
            self._term,
            self._trunc,
            {},
        )
