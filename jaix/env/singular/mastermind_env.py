import gymnasium as gym
from ttex.config import ConfigurableObject, Config
from jaix import LOGGER_NAME
import numpy as np
from typing import Tuple, Optional
import logging
from jaix.env.singular import SingularEnvironment

logger = logging.getLogger(LOGGER_NAME)


class MastermindEnvironmentConfig(Config):
    def __init__(
        self,
        num_slots_range: Tuple[int, int] = (10, 20),
        num_colours_range: Tuple[int, int] = (3, 5),
        sequential: bool = False,
        max_guesses: int = np.iinfo(np.int32).max,
    ):
        self.num_slots_range = num_slots_range
        self.num_colours_range = num_colours_range
        # Ensure correct rng format
        for rng in [num_colours_range, num_colours_range]:
            assert len(rng) == 2
            assert rng[0] <= rng[1]
        self.sequential = sequential
        self.max_guesses = max_guesses


class MastermindEnvironment(ConfigurableObject, SingularEnvironment):
    config_class = MastermindEnvironmentConfig

    def __init__(self, config: MastermindEnvironmentConfig, inst: int):
        ConfigurableObject.__init__(self, config)
        SingularEnvironment.__init__(self, inst)
        self._setup(config, inst)
        self.num_guesses = 0
        self.num_resets = 0

    def _setup(self, config: MastermindEnvironmentConfig, inst: int):
        np.random.seed(inst)
        self.num_slots = np.random.randint(
            low=config.num_slots_range[0], high=config.num_slots_range[1] + 1
        )
        self.num_colours = np.random.randint(
            low=config.num_colours_range[0],
            high=config.num_colours_range[1] + 1,
            size=self.num_slots,
        )
        self.action_space = gym.spaces.MultiDiscrete(self.num_colours)
        self.observation_space = gym.spaces.MultiDiscrete([self.num_slots + 1])

        self._solution = self.action_space.sample()

    def _get_info(self):
        return {
            "stop": self.stop(),
        }

    def stop(self):
        # TODO: potentially add fitness reached
        return self.num_guesses >= self.max_guesses

    def reset(
        self,
        *,
        seed: Optional[int] = None,
        options: Optional[dict] = None,
    ):
        """
        Resets the environment to an initial state,
        required before calling step.
        Returns the first agent observation for an episode and information,
        i.e. metrics, debug info.
        """
        self.num_resets += 1
        super().reset(seed=seed)
        return [self.num_slots], self._get_info()

    def step(self, x):
        """
        Updates an environment with actions returning the next agent observation,
        the reward for taking that actions,
        if the environment has terminated or truncated due to the latest action
        and information from the environment about the step,
        i.e. metrics, debug info.
        """
        self.num_guesses += 1
        x = np.asarray(x, dtype=self.action_space.dtype)
        # Obs is based on how many exact matches
        matches = x == self._solution
        if self.sequential:
            which_match = np.argwhere(matches == 0)
            if len(which_match) == 0:
                obs = self.num_slots
            else:
                obs = np.sum(matches[0 : which_match[0][0]])
        else:
            obs = np.sum(matches)
        # Minimisation
        obs = self.num_slots - obs
        terminated = obs == 0
        truncated = self.num_guesses >= self.max_guesses
        # observation, reward, terminated, truncated, info
        return [obs], obs, terminated, truncated, self._get_info()

    def render(self):
        logger.debug(self._get_info())

    def __str__(self):
        return f"MastermindEnvironment with {self.num_slots} slot, {self.num_colours} colours, solution {self._solution}"

    def close(self):
        pass
