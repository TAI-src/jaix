"""Factory to create environment from config and wrappers"""

from ttex.config import Config, ConfigurableObjectFactory as COF
from typing import Type, List, Tuple, Union, Dict
import gymnasium as gym
import logging

from jaix import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)


class WrappedEnvFactory:
    @staticmethod
    def wrap(
        env: gym.Env,
        wrappers: List[Tuple[Type[gym.Wrapper], Union[Config, Dict]]],
    ):
        wrapped_env = env
        for wrapper_class, wrapper_config in wrappers:
            logger.debug(f"Wrapping {env} with {wrapper_config} of {wrapper_class}")
            if isinstance(wrapper_config, Config):
                # Wrapper is a configurable object
                wrapped_env = COF.create(wrapper_class, wrapper_config, wrapped_env)
            else:
                # Assume config is a dict of keyword arguments
                wrapped_env = wrapper_class(wrapped_env, **wrapper_config)
        logger.debug(f"Wrapped env {wrapped_env}")
        return wrapped_env
