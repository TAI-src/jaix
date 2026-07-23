"""Factory to create environment from config and wrappers"""

import logging
from typing import cast

import gymnasium as gym
from ttex.config import (
    Config,
    ConfigurableObject,
)
from ttex.config import (
    ConfigFactory as CF,
)
from ttex.config import (
    ConfigurableObjectFactory as COF,
)

import jaix.utils.globals as globs

logger = logging.getLogger(globs.LOGGER_NAME)


class WrappedEnvFactory:
    @staticmethod
    def wrap(
        env: gym.Env,
        wrappers: list[tuple[type[gym.Wrapper], Config | dict]],
    ):
        wrapped_env = env
        for i, (wrapper_class, wrapper_config) in enumerate(wrappers):
            logger.debug(f"Wrapping {env} with {wrapper_config} of {wrapper_class}")
            if isinstance(wrapper_config, Config):
                # Wrapper is a configurable object and config is passed as object
                wrapped_env = COF.create(wrapper_class, wrapper_config, wrapped_env)
            elif issubclass(wrapper_class, ConfigurableObject):
                wrapper_class = cast(type[ConfigurableObject], wrapper_class)
                if (
                    len(wrapper_config) == 1
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
        logger.debug(f"Wrapped env {wrapped_env!s}")
        return wrapped_env
