"""Defines environment as in EC context"""

import gymnasium as gym
from gymnasium import spaces
import numpy as np
import logging

from ttex.config import ConfigurableObject, Config
from jaix.env.utils.problem.static_problem import StaticProblem
from typing import Optional
import jaix.utils.globals as globals
from jaix.env.singular.singular_environment import SingularEnvironment

logger = logging.getLogger(globals.LOGGER_NAME)


class ECEnvironmentConfig(Config):
    """EC Environment Config"""

    def __init__(self, budget_multiplier: int):
        """
        * budget_multiplier (int): Function budget = n*budget_multiplier
        """
        self.budget_multiplier = budget_multiplier


class ECEnvironment(ConfigurableObject, SingularEnvironment):
    """EC environment to run static problems with EC algorithms"""

    metadata = {"render_modes": ["ansi"], "render_fps": 4}
    config_class = ECEnvironmentConfig

    def __init__(
        self,
        config: ECEnvironmentConfig,
        func: StaticProblem,
        func_id: int = 0,  # This is just for string representation in SingularEnvironment
        inst: int = 0,  # This is just for string representation in SingularEnvironment
    ):
        ConfigurableObject.__init__(self, config)
        # TODO: get proper fun and inst id
        SingularEnvironment.__init__(self, func_id, inst)
        self.func = func
        # An action is a point in search space (x)
        self.action_space = spaces.Box(
            low=np.array(func.lower_bounds),
            high=np.array(func.upper_bounds),
            shape=(func.dimension,),
            dtype=np.float64,
        )
        # An observation are the objective values of the last action
        # f(x) = y
        self.observation_space = spaces.Box(
            low=np.array(func.min_values),
            high=np.array(func.max_values),
            shape=(func.num_objectives,),
            dtype=np.float64,
        )
        # Count how often the environment is reset
        # (corresponds to algorithms restarts +1 )
        self.num_resets = 0

    @property
    def suite_name(self) -> str:
        return f"ECEnvironment{type(self.func).__name__}"

    def _get_info(self):
        """Simple representation of environment state
        * num_resets: The number of time the environment has been reset
        * The number of evaluations left on the function
        """
        return {
            "num_resets": self.num_resets,
            "evals_left": self.func.evalsleft(self.budget_multiplier),
            "stop": self.stop(),
        }

    def stop(self):
        return self.func.stop(self.budget_multiplier)

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
            raise ValueError("EC environments are always online")
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
        raw_fitness, clean_fitness = self.func(x)
        terminated = self.func.final_target_hit()
        truncated = (
            self.func.evalsleft(self.budget_multiplier) <= 0
        )  # or not self.action_space.contains(x)
        # Not truncating if out of bounds for EC
        # If algorithm wants to implement that, that is fine
        # TODO: potential constraint handling strategies
        # if out out bounds

        # observation, reward, terminated, truncated, info
        # TODO: reward cannot be multi-dimensional
        return raw_fitness, None, terminated, truncated, self._get_info()

    def render(self):
        """
        Renders the environments to help visualise what the agent see,
        examples modes are “human”, “rgb_array”, “ansi” for text.
        """
        logger.debug(self._get_info())

    def close(self):
        """
        Closes the environment, important when external software is used,
        i.e. pygame for rendering, databases
        """
        return self.func.close()

    @property
    def name(self):
        return f"ECEnvironment/{self.func.__class__.__name__}"
