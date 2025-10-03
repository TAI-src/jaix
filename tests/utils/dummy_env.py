import gymnasium as gym
from gymnasium import spaces
from typing import Optional
from gymnasium.utils.env_checker import check_env
from ttex.config import ConfigurableObject, Config
from jaix.env.singular.singular_environment import SingularEnvironment


class DummyEnvConfig(Config):
    def __init__(
        self,
        dimension: int = 3,
        num_objectives: int = 1,
        action_space=None,
        observation_space=None,
        reward_space=None,
    ):
        self.dimension = dimension
        self.num_objectives = num_objectives
        self.action_space = action_space
        self.observation_space = observation_space
        self.reward_space = reward_space


class DummyConfEnv(ConfigurableObject, SingularEnvironment):
    config_class = DummyEnvConfig

    @staticmethod
    def info(*args, **kwargs):
        return {
            "funcs": [0, 1],
            "insts": list(range(15)),
            "val": 13,
        }

    def __init__(self, config: DummyEnvConfig, func: int = 0, inst: int = 1):
        ConfigurableObject.__init__(self, config)
        SingularEnvironment.__init__(self, func, inst)
        if self.action_space is None:
            self.action_space = spaces.Box(low=-5, high=5, shape=(self.dimension,))
        if self.observation_space is None:
            self.observation_space = spaces.Box(
                low=0, high=100, shape=(self.num_objectives,)
            )
        if self.reward_space is None:
            self.reward_space = spaces.Box(low=0, high=5)
        self._trunc = False
        self._term = False
        self._stop = False
        self.inst = inst
        self.func = func
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
        return self.observation_space.sample(), DummyEnv.info()

    def step(self, x):
        return (
            self.observation_space.sample(),
            self.reward_space.sample()[0],
            self._term,
            self._trunc,
            DummyEnv.info(),
        )

    def stop(self):
        return self._stop


class DummyEnv(DummyConfEnv):
    def __init__(
        self,
        dimension=3,
        num_objectives=1,
        action_space=None,
        observation_space=None,
        reward_space=None,
    ):
        config = DummyEnvConfig(
            dimension=dimension,
            num_objectives=num_objectives,
            action_space=action_space,
            observation_space=observation_space,
            reward_space=reward_space,
        )
        DummyConfEnv.__init__(self, config)


def test_dummy_env():
    check_env(DummyEnv(), skip_render_check=True)
