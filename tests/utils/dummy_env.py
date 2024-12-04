import gymnasium as gym
from gymnasium import spaces
from typing import Optional
from gymnasium.utils.env_checker import check_env


class DummyEnv(gym.Env):
    def __init__(self, dimension=3, num_objectives=1):
        self.action_space = spaces.Box(low=-5, high=5, shape=(dimension,))
        self.observation_space = spaces.Box(low=0, high=100, shape=(num_objectives,))
        self.reward_space = spaces.Box(low=0, high=5)
        self._trunc = False
        self._term = False
        self._stop = False
        super().__init__()

    def reset(
        self,
        *,
        seed: Optional[int] = None,
        options: Optional[dict] = None,
    ):
        super().reset(seed=seed)
        self.observation_space.seed(seed)
        self.action_space.seed(seed)
        self.reward_space.seed(seed)
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

    def stop(self):
        return self._stop

    def __str__(self):
        return "DummyEnv"


def test_dummy_env():
    check_env(DummyEnv(), skip_render_check=True)
