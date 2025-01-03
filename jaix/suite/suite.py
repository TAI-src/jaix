from enum import Enum
from typing import Optional, Type
from ttex.config import ConfigurableObject, ConfigurableObjectFactory as COF, Config
import gymnasium as gym
import random as rnd
import logging
from jaix import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)


class AggType(Enum):
    NONE = 0
    INST = 1


class SuiteConfig(Config):
    def __init__(
        self,
        env_class: Type[gym.Env],
        env_config: Config,
        num_instances: int,
        num_agg_instances: int,
    ):
        self.env_class = env_class
        self.env_config = env_config
        self.num_instances = num_instances
        self.num_agg_instances = num_agg_instances


class Suite(ConfigurableObject):
    config_class = SuiteConfig

    def _get_env(self, inst):
        return COF.create(self.env_class, self.env_config, inst)

    def get_envs(self):
        for inst in range(self.num_instances):
            env = self._get_env(inst)
            yield env

    def get_agg_envs(self, agg_type: AggType, seed: Optional[int] = None):
        logger.debug(f"Getting environments with seed {seed}")
        seeds = rnd.Random(seed).sample(range(100, 999999), self.num_instances)
        for i_seed in seeds:
            logger.debug(f"Shuffle seed is {i_seed}")
            def_envs = [self._get_env(inst) for inst in range(self.num_agg_instances)]
            logger.debug(f"Initialised {len(def_envs)} sub-environments")
            # shuffle instance ids
            inst_ids = rnd.Random(i_seed).sample(
                range(self.num_agg_instances), self.num_agg_instances
            )
            envs = [def_envs[i] for i in inst_ids]
            logger.debug(f"Returning {envs}")
            yield envs
