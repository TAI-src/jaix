import gymnasium as gym
from ttex.config import ConfigurableObject, Config
from jaix import LOGGER_NAME
import numpy as np
from typing import Tuple, Optional
import logging

logger = logging.getLogger(LOGGER_NAME)


class MastermindEnvironmentConfig(Config):
    def __init__(
        self,
        num_slots_range: Tuple[int, int],
        num_colours_range: Tuple[int, int],
        sequential: bool,
        num_guesses: int = np.inf,
    ):
        self.num_slots_range = num_slots_range
        self.num_colours_range = num_colours_range
        # Ensure correct rng format
        for rng in [num_colours_range, num_colours_range]:
            assert len(rng) == 2
            assert rng[0] <= rng[1]
        self.sequential = sequential
        self.num_guesses = num_guesses


class MastermindEnvironment(ConfigurableObject, gym.Env):
    config_class = MastermindEnvironmentConfig

    def __init__(self, config: MastermindEnvironmentConfig, inst: int):
        ConfigurableObject.__init__(self, config)
        self._setup(config, inst)
        self.action_space = gym.spaces.MultiDiscrete(
            [self.num_colours] * self.num_guesses
        )
        self.observation_space = gym.spaces.Discrete(self.num_slots)
        self.guesses_done = 0

    def _setup(self, config: MastermindEnvironmentConfig, inst: int):
        np.random.seed(inst)
        self.num_colours = np.random.randint(
            low=config.num_colours_range[0], high=config.num_colours_range[1]
        )
        self.num_slots = np.random.randint(
            low=config.num_slots_range[0], high=config.num_slots_range[1]
        )
        self.solution = np.random.randint(self.num_colours, size=self.num_slots)

    def _get_info(self):
        return {"stop": self.stop(), "guesses_done": self.guesses_done}

    def stop(self):
        # TODO: potentially add fitness reached
        return self.guesses_done >= self.num_guesses

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
        if options is None or "online" not in options or not options["online"]:
            # We only do partial resets for ec, so still "online"
            raise ValueError("Mastermind environments are always online")
        self.num_resets += 1
        return None, self._get_info()

    def step(self, x):
        """
        Updates an environment with actions returning the next agent observation,
        the reward for taking that actions,
        if the environment has terminated or truncated due to the latest action
        and information from the environment about the step,
        i.e. metrics, debug info.
        """
        x = np.asarray(x, dtype=self.action_space.dtype)
        # Obs is based on how many exact matches
        matches = x == self.solution
        if self.sequential:
            which_match = np.argwhere(matches)
            if len(which_match) == 0:
                obs = 0
            else:
                obs = np.sum(matches[0 : which_match[0]])
        else:
            obs = np.sum(matches)
        terminated = obs == self.num_slots
        truncated = self.guesses_done >= self.num_guesses
        # observation, reward, terminated, truncated, info
        return obs, obs, terminated, truncated, self._get_info()

    def render(self):
        logger.debug(self._get_info())

    def __str__(self):
        return f"MastermindEnvironment with {self.num_slots} slot, {self.num_colours} colours"

    def close(self):
        pass
