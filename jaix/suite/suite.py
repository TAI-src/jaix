from enum import Enum
from typing import Optional, Type, List, Union, Tuple
from ttex.config import ConfigurableObject, ConfigurableObjectFactory as COF, Config
import gymnasium as gym
import random as rnd
import logging
from jaix import LOGGER_NAME
import itertools


logger = logging.getLogger(LOGGER_NAME)


class AggType(Enum):
    NONE = 0
    INST = 1


class SuiteConfig(Config):
    def __init__(
        self,
        env_class: Type[gym.Env],
        env_config: Config,
        functions: Optional[List[int]] = None,
        instances: Optional[List[int]] = None,
        agg_instances: Optional[Union[List[int], List[Tuple[int]], int]] = None,
    ):
        self.env_class = env_class
        self.env_config = env_config
        env_info = {}
        if hasattr(env_class, "info"):
            env_info = env_class.info(env_config)
        # TODO: Better error messages if no functions / instances passed
        # and info does not exist
        self.functions = env_info["funcs"] if functions is None else functions
        self.instances = env_info["insts"] if instances is None else instances
        instance_permutations = list(itertools.permutations(self.instances))
        if agg_instances is None:
            self.agg_instances = instance_permutations
        elif isinstance(agg_instances, int):
            assert agg_instances < len(instance_permutations)
            self.agg_instances = [instance_permutations[agg_instances]]
        elif isinstance(agg_instances, list) and all(
            [isinstance(i, int) for i in agg_instances]
        ):
            self.agg_instances = [instance_permutations[i] for i in agg_instances]  # type: ignore
        else:
            self.agg_instances = agg_instances  # type: ignore


class Suite(ConfigurableObject):
    config_class = SuiteConfig

    def _get_env(self, func, inst):
        return COF.create(self.env_class, self.env_config, func, inst)

    def get_envs(self):
        for func in self.functions:
            for inst in self.instances:
                logger.warning(
                    f"Getting environment for function {func} and instance {inst}"
                )
                env = self._get_env(func, inst)
                yield env

    def get_agg_envs(self, agg_type: AggType, seed: Optional[int] = None):
        logger.debug(f"Getting environments with seed {seed}")
        if agg_type != AggType.INST:
            raise NotImplementedError("Only INST aggregation is supported")
        for func in self.functions:
            for agg_inst in self.agg_instances:
                envs = [self._get_env(func, inst) for inst in agg_inst]
                logger.debug(f"Returning {envs}")
                yield envs
