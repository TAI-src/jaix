from ttex.config import ConfigurableObject, Config
import gymnasium as gym
import numpy as np
from tabrepo.repository.evaluation_repository import (
    load_repository,
)
from jaix.env.utils.hpo.tabrepo_adapter import TaskType, TabrepoAdapter
from typing import Optional, List, Tuple, Dict
import jaix.utils.globals as globals
from jaix.env.singular.singular_environment import SingularEnvironment
from collections import defaultdict
import json

# TODO: Introduce ensembles at some point
import logging

logger = logging.getLogger(globals.LOGGER_NAME)


class HPOEnvironmentConfig(Config):
    def __init__(
        self,
        task_type: TaskType,
        training_budget: int = np.iinfo(np.int32).max,
        repo_name: str = "D244_F3_C1530_30",
        cache: bool = True,
        target_rank: int = 1,
    ):
        self.training_budget = training_budget
        self.repo = load_repository(repo_name, load_predictions=True, cache=cache)
        self.task_type = task_type
        self.target_rank = target_rank


class HPOEnvironment(ConfigurableObject, SingularEnvironment):
    config_class = HPOEnvironmentConfig

    @staticmethod
    def info(config: HPOEnvironmentConfig):
        datasets = TabrepoAdapter.get_dataset_names(config.repo, config.task_type)
        num_funcs = len(datasets)
        num_insts = 3  # TODO: This is an assumption
        return {"funcs": list(range(num_funcs)), "insts": list(range(num_insts))}

    def __init__(
        self,
        config: HPOEnvironmentConfig,
        func: int,
        inst: int,
    ):
        ConfigurableObject.__init__(self, config)
        SingularEnvironment.__init__(self, func, inst)
        self.tabrepo_adapter = TabrepoAdapter(
            self.repo, self.task_type, dataset_idx=func, fold=inst
        )
        # An action is the index of a config
        # which is basically the type of model chosen
        # for the ensemble
        # TODO: proper config space with actual hyperparameters
        self.n = len(self.tabrepo_adapter.configs)
        self.action_space = gym.spaces.MultiBinary(self.n)
        # Observation is the validation error of the last config
        self.observation_space = gym.spaces.Box(
            low=0,
            high=self.tabrepo_adapter.max_rank,
            shape=(1,),
            dtype=np.float64,
        )
        self.training_time = 0
        self.num_resets = 0
        self.ensembles = defaultdict(
            list
        )  # type: Dict[float, List[Tuple[List[int], float]]]

    def _get_info(self):
        # TODO: don't send full ensembles mid-way
        return {
            "dataset": self.tabrepo_adapter.metadata,
            "stop": self.stop(),
            "env_step": self.training_time,
            "ensembles": json.dumps(self.ensembles),
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
        config_ids = np.where(x)[0]
        obs, time_train_s = self.tabrepo_adapter.evaluate_ensemble(config_ids)
        logger.debug(
            f"Action {config_ids} resulted in obs {obs} with time {time_train_s}"
        )
        # Record
        self.ensembles[obs].append((config_ids.tolist(), time_train_s))

        self.training_time += time_train_s
        terminated = obs < self.target_rank
        truncated = self.stop()
        info = self._get_info()
        info["time_train_s"] = time_train_s
        return [obs], None, terminated, truncated, info

    def render(self):
        """
        Renders the environments to help visualise what the agent see,
        examples modes are “human”, “rgb_array”, “ansi” for text.
        """
        logger.debug(self._get_info())
