import gymnasium as gym
from gymnasium import spaces
from typing import Optional
from gymnasium.utils.env_checker import check_env
from ttex.config import ConfigurableObject, Config


class DummyEnvConfig(Config):
    def __init__(
        self,
        dimension: int = 3,
        num_objectives: int = 1,
    ):
        self.dimension = dimension
        self.num_objectives = num_objectives


class DummyConfEnv(gym.Env, ConfigurableObject):
    config_class = DummyEnvConfig

    def __init__(self, config: DummyEnvConfig, inst: int = 1):
        ConfigurableObject.__init__(self, config)
        self.action_space = spaces.Box(low=-5, high=5, shape=(self.dimension,))
        self.observation_space = spaces.Box(
            low=0, high=100, shape=(self.num_objectives,)
        )
        self.reward_space = spaces.Box(low=0, high=5)
        self._trunc = False
        self._term = False
        self._stop = False
        self.inst = inst
        self.id = 42

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


class DummyEnv(DummyConfEnv):
    def __init__(self, dimension=3, num_objectives=1):
        config = DummyEnvConfig(dimension=dimension, num_objectives=num_objectives)
        DummyConfEnv.__init__(self, config)


def test_dummy_env():
    check_env(DummyEnv(), skip_render_check=True)
