from enum import Enum
from typing import Optional, Type, List, Union, Tuple, cast
from ttex.config import ConfigurableObject, ConfigurableObjectFactory as COF, Config
import gymnasium as gym
import random as rnd
import logging
from jaix.utils.globals import LOGGER_NAME
import itertools
import math
import numpy as np

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
        comp_env_num: Optional[int] = None,
        agg_instances: Optional[Union[List[int], List[Tuple[int, ...]], int]] = None,
        seed: Optional[int] = None,
    ):
        logger.debug("Creating SuiteConfig")
        self.env_class = env_class
        self.env_config = env_config
        env_info = {}
        if hasattr(env_class, "info"):
            env_info = env_class.info(env_config)
        # TODO: Better error messages if no functions / instances passed
        # and info does not exist
        self.functions = env_info["funcs"] if functions is None else functions
        self.instances = env_info["insts"] if instances is None else instances
        # TODO: Allow function / instance combinations
        # generate instance permuations of length comp_env_num
        comp_env_num = len(self.instances) if comp_env_num is None else comp_env_num
        instance_permutations = itertools.permutations(self.instances, comp_env_num)

        self.agg_instances: List[Tuple[int, ...]]
        if agg_instances is None:
            # Nothing passed, use all permutations of instances
            if len(self.instances) > 5:
                logger.warning(
                    "No aggregation instances passed, using all permutations of instances. "
                    "This may take a long time for large instance sets or crash the system."
                )
            self.agg_instances = list(instance_permutations)
        elif isinstance(agg_instances, int):
            # Integer n passed, use n random permutations
            # TODO: Make this seedable
            assert agg_instances > 0, "agg_instances must be a positive integer"
            if len(self.instances) < 5:
                logger.warning(
                    "Using random permutations of instances. Inefficient for small instance sets."
                )
            inst_tuples: List[Tuple[int, ...]] = [
                tuple(
                    np.random.choice(self.instances, size=comp_env_num, replace=False)
                )
                for _ in range(agg_instances)
            ]
            while len(set(inst_tuples)) < agg_instances:
                # Ensure we have unique tuples
                inst_tuples.append(
                    tuple(
                        np.random.choice(
                            self.instances, size=comp_env_num, replace=False
                        )
                    )
                )
            # expecting normal ints, so reformatting
            self.agg_instances = [
                tuple([int(x) for x in perm]) for perm in set(inst_tuples)
            ]
        elif isinstance(agg_instances, list):
            if all([isinstance(i, int) for i in agg_instances]):
                agg_instances = cast(
                    List[int], agg_instances
                )  # otherwise mypy complains
                assert all(
                    i >= 0 for i in agg_instances
                ), "agg_instances must be a list of non-negative integers"
                assert all(
                    i < math.factorial(len(self.instances)) for i in agg_instances
                ), "agg_instances must be a list of integers less than the number of permutations (0 to n!)"
                # List of integers passed (the desired indices). iterate to find the permutations
                if len(self.instances) > 5:
                    logger.warning(
                        "Filtering out indices. This may take a long time for large instance sets"
                    )
                self.agg_instances = []
                for i, perm in enumerate(instance_permutations):
                    if i in agg_instances:
                        self.agg_instances.append(perm)
            elif all([isinstance(i, tuple) for i in agg_instances]):
                agg_instances = cast(
                    List[Tuple[int, ...]], agg_instances
                )  # otherwise mypy complains
                assert all(
                    min(i) >= 0 for i in agg_instances
                ), "agg_instances must be a list of tuples with non-negative integers"
                assert all(
                    max(i) < len(self.instances) for i in agg_instances
                ), "agg_instances must be a list of tuples with integers less than the number of instances"
                # List of tuples passed (the desired permutations). Use those directly
                self.agg_instances = agg_instances
        else:
            # Invalid type passed, raise error
            raise ValueError(
                "agg_instances must be None, an integer, or a list of tuples"
            )
        logger.debug(f"SuiteConfig created with {self.__dict__}")


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
