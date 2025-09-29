from ttex.config import (
    ConfigurableObject,
    ConfigurableObjectFactory as COF,
    Config,
)
from jaix.env.utils.switching_pattern.switching_pattern import (
    SwitchingPattern,
)
from jaix.env.wrapper.auto_reset_wrapper import (
    AutoResetWrapper,
    AutoResetWrapperConfig,
)
from jaix.env.wrapper.wrapped_env_factory import WrappedEnvFactory as WEF
from typing import Dict, Type, List, Optional, Callable, Any, TypeVar, Tuple, Union
import gymnasium as gym
from gymnasium import spaces
from functools import wraps
from jaix.env.composite.composite_environment import CompositeEnvironment

import logging
from jaix.utils.globals import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
FuncT = TypeVar("FuncT", bound=Callable[..., Any])


class SwitchingEnvironmentConfig(Config):
    def __init__(
        self,
        switching_pattern_class: Type[SwitchingPattern],
        switching_pattern_config: Config,
        real_time: bool,
        auto_reset_wrapper_config: Optional[AutoResetWrapperConfig] = None,
    ):
        self.switching_pattern_class = switching_pattern_class
        self.switching_pattern_config = switching_pattern_config
        self.real_time = real_time
        if auto_reset_wrapper_config is None:
            self.auto_reset_wrapper_config = AutoResetWrapperConfig()
        else:
            self.auto_reset_wrapper_config = auto_reset_wrapper_config


class SwitchingEnvironment(ConfigurableObject, CompositeEnvironment):
    """Environment that dynamically switches between a list of environments depending on time"""

    config_class = SwitchingEnvironmentConfig

    def __init__(
        self,
        config: SwitchingEnvironmentConfig,
        env_list: List[gym.Env],
    ):
        ConfigurableObject.__init__(self, config)
        wrappers = [
            (AutoResetWrapper, self.auto_reset_wrapper_config),
        ]  # type: List[Tuple[Type[gym.Wrapper], Union[Config, Dict]]]
        self.env_list = [WEF.wrap(env, wrappers) for env in env_list]
        self._current_env = 0

        self.pattern_switcher = COF.create(
            self.switching_pattern_class,
            self.switching_pattern_config,
            num_choices=len(self.env_list),
        )
        if self.real_time:
            # save current wallclocktime etc
            raise NotImplementedError()
        else:
            self._timer = 0
        self.steps_counter = 0
        self._stopped = False
        self._set_spaces()
        self.constant_dim = CompositeEnvironment.const_dim(env_list)

    def _set_spaces(self):
        env = self.env_list[self._current_env]
        self.action_space = env.action_space
        switcher_space = spaces.Discrete(len(self.env_list))
        self.observation_space = spaces.Tuple((switcher_space, env.observation_space))
        logger.debug(f"Action space: {self.action_space}")
        logger.debug(f"Observation space: {self.observation_space}")

    def update_env(func: FuncT):
        @wraps(func)
        def decorator_func(self, *args: Any, **kwargs: Any) -> Any:
            self._update_current_env()
            return func(self, *args, **kwargs)

        return decorator_func

    def _get_info(self):
        # Not updating here to avoid double-update.
        # Function private
        return {
            "steps_counter": self.steps_counter,
            "timer": self._timer,
            "current_env": self._current_env,
        }

    @update_env
    def step(self, action):
        obs, r, term, trunc, info = self.env_list[self._current_env].step(action)
        obs = (self._current_env, obs)

        if not self.real_time:
            # Counting next function evaluation
            self._timer = self._timer + 1

        self.steps_counter = self.steps_counter + 1

        info["meta"] = self._get_info()
        logger.debug(info)

        return obs, r, term, self._stopped, info

    def reset(
        self,
        *,
        seed: Optional[int] = None,
        options: Optional[dict] = None,
    ):
        if options is None:
            options = {}
        if "online" in options and options["online"]:
            env_obs, info = self.env_list[self._current_env].reset(
                seed=seed, options=options
            )
            logger.debug("Online reset")
        else:
            self._timer = 0
            self._current_env = 0
            self.pattern_switcher.reset(seed=seed)
            # Reset all environments in the list
            for env in reversed(self.env_list):
                # Reversed order so return value is of current environment
                # which is the first one in the list (after reset)
                env_obs, info = env.reset(seed=seed, options=options)
            logger.debug("Full reset")
        self._update_current_env()
        obs = (self._current_env, env_obs)
        info["meta"] = self._get_info()
        return obs, info

    @update_env
    def render(self):
        return self.env_list[self._current_env].render()

    def _update_current_env(self):
        # Increase real-time timer
        if self.real_time:
            # TODO: get time now
            raise NotImplementedError()

        # Get current environment and return accordingly
        valid_envs = [not env.stop() for env in self.env_list]
        new_env = self.pattern_switcher.switch(self._timer, valid=valid_envs)
        # only update env if not invalid (-1)
        self._stopped = True if new_env < 0 else False
        updated = False
        new_env = max(0, new_env)
        if new_env != self._current_env:
            updated = True
        self._current_env = new_env
        if updated:
            self._set_spaces()

    def _stop(self):
        stopped_envs = [env.stop() for env in self.env_list]
        return all(stopped_envs) or self._stopped

    @update_env
    def stop(self):
        return self._stopped

    def close(self):
        return [env.close() for env in self.env_list]

    @update_env
    def __getattr__(self, name):
        # Pass through anything else that is not overriden
        # This function is only called after it is not found in self
        return getattr(self.env_list[self._current_env].unwrapped, name)
