import gymnasium as gym
from gymnasium import spaces
from typing import Optional


class DummyEnv(gym.Env):
    def __init__(self):
        self.action_space = spaces.Discrete(2)
        self.observation_space = spaces.Discrete(5)
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
