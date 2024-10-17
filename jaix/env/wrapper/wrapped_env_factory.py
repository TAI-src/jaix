"""Factory to create environment from config and wrappers """
from ttex.config import Config, ConfigurableObjectFactory as COF
from typing import Type, List, Tuple, Union, Dict
import gymnasium as gym


class WrappedEnvFactory:
    @staticmethod
    def wrap(
        env: gym.Env,
        wrappers: List[Tuple[Type[gym.Wrapper], Union[Config, Dict]]],
    ):
        wrapped_env = env
        for wrapper_class, wrapper_config in wrappers:
            if isinstance(wrapper_config, Config):
                # Wrapper is a configurable object
                wrapped_env = COF.create(wrapper_class, wrapper_config, wrapped_env)
            else:
                # Assume config is a dict of keyword arguments
                wrapped_env = wrapper_class(wrapped_env, **wrapper_config)
        return wrapped_env
