"""Factory to create environment from config and wrappers"""

from ttex.config import (
    Config,
    ConfigurableObjectFactory as COF,
    ConfigurableObject,
    ConfigFactory as CF,
)
from typing import Type, List, Tuple, Union, Dict
import gymnasium as gym
import logging

from jaix.utils.globals import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)


class WrappedEnvFactory:
    @staticmethod
    def wrap(
        env: gym.Env,
        wrappers: List[Tuple[Type[gym.Wrapper], Union[Config, Dict]]],
    ):
        wrapped_env = env
        for i, (wrapper_class, wrapper_config) in enumerate(wrappers):
            logger.debug(f"Wrapping {env} with {wrapper_config} of {wrapper_class}")
            if isinstance(wrapper_config, Config):
                # Wrapper is a configurable object and config is passed as object
                wrapped_env = COF.create(wrapper_class, wrapper_config, wrapped_env)
            elif (
                issubclass(wrapper_class, ConfigurableObject)
                and len(wrapper_config) == 1
                and str(
                    f"{wrapper_class.config_class.__module__}.{wrapper_class.config_class.__qualname__}"
                )
                in wrapper_config
            ):
                # Wrapper is a configurable object and config is passed as dict
                config_object = CF.from_dict(wrapper_config, context=globals())
                wrapped_env = COF.create(wrapper_class, config_object, wrapped_env)
                wrappers[i] = (
                    wrapper_class,
                    config_object,
                )  # Update to use config object
            else:
                # Assume config is a dict of keyword arguments
                wrapped_env = wrapper_class(wrapped_env, **wrapper_config)
        logger.debug(f"Wrapped env {wrapped_env}")  #
        return wrapped_env
