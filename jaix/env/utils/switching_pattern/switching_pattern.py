from ttex.config import ConfigurableObject, ConfigurableObjectFactory as COF, Config
import math
from abc import abstractmethod
from gymnasium import spaces
from typing import List, Optional


class SwitchingPattern(ConfigurableObject):
    def __init__(self, config: Config, num_choices: int):
        ConfigurableObject.__init__(self, config)
        self.num_choices = num_choices

    @abstractmethod
    def reset(self, seed=None):
        # reset so the surrounding environment can reset
        raise NotImplementedError()

    @abstractmethod
    def switch(self, t: float, valid: List[bool]) -> int:
        # Return which env index we are switching to
        raise NotImplementedError()


class SeqRegSwitchingPatternConfig(Config):
    def __init__(
        self,
        wait_period: float,
        carry_over: bool = False,
    ):
        self.wait_period = wait_period
        self.carry_over = carry_over
        self.offset = 0


class SeqRegSwitchingPattern(SwitchingPattern):
    config_class = SeqRegSwitchingPatternConfig

    def reset(self, seed=None):
        self.offset = 0

    def switch(self, t: float, valid: Optional[List[bool]] = None) -> int:
        assert valid is None or self.num_choices == len(valid)
        valid = [True] * self.num_choices if valid is None else valid
        env_num = math.floor((t + self.offset) / self.wait_period)
        if env_num >= self.num_choices:
            # Out of environments, return -1
            return -1
        if not valid[env_num]:
            # Chosen environment stopped, get the next in line
            # If none is left, return -1
            env_num = next(
                (env for env in range(env_num + 1, self.num_choices) if valid[env]), -1
            )
            if not self.carry_over:
                current = env_num if env_num >= 0 else self.num_choices
                # Add an internal offset as remaining budget cannot be carried over
                self.offset = current * self.wait_period - t
                # TODO: better tests for offset
        return env_num


class SeqForcedSwitchingPatternConfig(Config):
    def __init__(
        self,
    ):
        self.current = 0


class SeqForcedSwitchingPattern(SwitchingPattern):
    config_class = SeqForcedSwitchingPatternConfig

    def reset(self, seed=None):
        self.current = 0

    def switch(self, t: float, valid: List[bool]) -> int:
        if not valid[self.current]:
            self.current = next(
                (
                    env
                    for env in range(self.current + 1, self.num_choices)
                    if valid[env]
                ),
                -1,
            )
        if self.current >= self.num_choices:
            return -1
        return self.current
