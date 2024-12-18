from ttex.config import ConfigurableObject, Config
import gymnasium as gym
import numpy as np
from tabrepo.repository.evaluation_repository import EvaluationRepository
from jaix.env.utils.hpo import TaskType, TabrepoAdapter
from typing import Optional
from jaix import LOGGER_NAME

# TODO: Introduce ensembles at some point
import logging

logger = logging.getLogger(LOGGER_NAME)


class HPOEnvironmentConfig(Config):
    def __init__(self, training_budget: int):
        self.training_budget = training_budget


class HPOEnvironment(ConfigurableObject, gym.Env):
    config_class = HPOEnvironmentConfig

    def __init__(
        self,
        config: HPOEnvironmentConfig,
        repo: EvaluationRepository,
        task_type: TaskType,
        inst: int,
    ):
        ConfigurableObject.__init__(self, config)
        self.tabrepo_adapter = TabrepoAdapter(repo, task_type, inst)
        # An action is the index of a config
        # TODO: proper config space with actual hyperparameters
        self.action_space = gym.spaces.Discrete(len(self.tabrepo_adapter.configs))
        # Observation is the validation error of the last config
        self.observation_space = gym.spaces.Box(
            low=np.array([self.tabrepo_adapter.metadata["min_error_val"]]),
            high=np.array([self.tabrepo_adapter.metadata["max_error_val"]]),
            shape=(1,),
            dtype=np.float64,
        )
        self.training_time = 0
        self.num_resets = 0

    def _get_info(self):
        return {
            "dataset": self.tabrepo_adapter.metadata,
            "training_time": self.training_time,
            "stop": self.stop(),
        }

    def stop(self):
        return self.training_time >= self.training_budget

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
            raise ValueError("HPO environments are always online")
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
        obs, time_train_s = self.tabrepo_adapter.evaluate(x)
        self.training_time += time_train_s
        terminated = False
        truncated = self.stop()
        return [obs], obs, terminated, truncated, self._get_info()

    def render(self):
        """
        Renders the environments to help visualise what the agent see,
        examples modes are “human”, “rgb_array”, “ansi” for text.
        """
        logger.debug(self._get_info())

    def __str__(self):
        return f"HPO Environment {str(self.tabrepo_adapter)}"
