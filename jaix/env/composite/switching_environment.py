from ttex.config import (
    ConfigurableObject,
    ConfigurableObjectFactory as COF,
    Config,
)
from jaix.env.utils.switching_pattern import (
    SwitchingPattern,
)
from jaix.env.wrapper import AutoResetWrapper, AnyFitWrapper
from typing import Type, List, Optional
import gymnasium as gym
from gymnasium import spaces
import numpy as np

import logging

logger = logging.getLogger("DefaultLogger")


class SwitchingEnvironmentConfig(Config):
    def __init__(
        self,
        switching_pattern_class: Type[SwitchingPattern],
        switching_pattern_config: Config,
        real_time: bool,
        next_after_resets: int = np.inf,
    ):
        self.switching_pattern_class = switching_pattern_class
        self.switching_pattern_config = switching_pattern_config
        self.real_time = real_time
        self.next_after_resets = next_after_resets
        assert self.next_after_resets > 0


class SwitchingEnvironment(ConfigurableObject, gym.Env):
    """Environment that dynamically switches between a list of environments depending on time"""

    config_class = SwitchingEnvironmentConfig

    def __init__(
        self,
        config: SwitchingEnvironmentConfig,
        env_list: List[gym.Env],
        observation_space: spaces.Space,
        action_space: spaces.Space,
    ):
        ConfigurableObject.__init__(self, config)
        self.env_list = [AutoResetWrapper(AnyFitWrapper(env)) for env in env_list]
        self._current_env = 0

        self.pattern_switcher = COF.create(
            self.switching_pattern_class, self.switching_pattern_config
        )
        if self.real_time:
            # save current wallclocktime etc
            raise NotImplementedError()
        else:
            self._timer = 0
        switcher_space = spaces.Discrete(self.switching_pattern_config.num_choices)
        self.steps_counter = 0

        self.action_space = action_space
        self.observation_space = spaces.Tuple((switcher_space, observation_space))
        self._stopped = False

    def update_env(func):
        def decorator_func(self, *args, **kwargs):
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

        return obs, r, term, trunc or self._stop(), info

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
        else:
            self._timer = 0
            self._current_env = 0
            self.pattern_switcher.reset(seed=seed)
            # Reset all environments in the list
            for env in reversed(self.env_list):
                # Reversed order so return value is of current environment
                # which is the first one in the list (after reset)
                env_obs, info = env.reset(seed=seed, options=options)
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
        valid_envs = [not env.stop(self.next_after_resets) for env in self.env_list]
        new_env = self.pattern_switcher.switch(self._timer, valid=valid_envs)
        # only update env if not invalid (-1)
        self._stopped = True if new_env < 0 else False
        self._current_env = max(0, new_env)

    def _stop(self):
        stopped_envs = [env.stop(self.next_after_resets) for env in self.env_list]
        return all(stopped_envs) or self._stopped

    @update_env
    def stop(self):
        return self._stop()

    def close(self):
        rec_files = [env.close() for env in self.env_list]
        self.closed = True
        return rec_files

    @update_env
    def __getattr__(self, name):
        # Pass through anything else that is not overriden
        # This function is only called after it is not found in self
        return getattr(self.env_list[self._current_env].unwrapped, name)