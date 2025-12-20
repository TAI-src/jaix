import gymnasium as gym
import numpy as np
import math
from jaix.env.wrapper.value_track_wrapper import (
    ValueTrackWrapper,
)
from ttex.config import ConfigurableObject, Config
from typing import Optional, Any, Dict, Tuple
from enum import Enum


class ImprovementType(Enum):
    OVER_FIRST = 1
    OVER_BEST = 2
    OVER_LAST = 3
    BEST_SINCE_FIRST = 4


class ImprovementRewardWrapperConfig(Config):
    def __init__(
        self,
        state_eval: str = "obs0",
        is_min: bool = True,
        transform: bool = True,
        passthrough: bool = True,
        imp_type=ImprovementType.BEST_SINCE_FIRST,
    ):
        self.state_eval = state_eval
        self.is_min = is_min
        self.transform = transform
        self.passthrough = passthrough
        self.imp_type = imp_type


class ImprovementRewardWrapper(ConfigurableObject, ValueTrackWrapper):
    config_class = ImprovementRewardWrapperConfig
    """
    A reward wrapper that transforms the state evaluation into a reward based on improvement.
    """

    def __init__(self, config: ImprovementRewardWrapperConfig, env: gym.Env):
        ConfigurableObject.__init__(self, config)
        ValueTrackWrapper.__init__(
            self, env, config.state_eval, config.is_min, config.passthrough
        )

    @staticmethod
    def _pos_log_scale_axis(val: float) -> float:
        if abs(val) <= 1:
            return 0
        ret = math.copysign(1, val) * np.log10(abs(val))
        return ret

    def _compute_imp(self, comp_val: float, val: float) -> float:
        if self.is_min:
            imp = comp_val - val
        else:
            imp = val - comp_val
        return max(0, imp)

    def _get_improvement(self, val: float) -> float:
        # On first step, just return 0 improvement
        if self.first_val is None or self.best_val is None or self.last_val is None:
            return 0.0
        if self.imp_type == ImprovementType.OVER_FIRST:
            comp_val = self.first_val
        elif self.imp_type == ImprovementType.OVER_BEST:
            comp_val = self.best_val
        elif self.imp_type == ImprovementType.OVER_LAST:
            comp_val = self.last_val
        elif self.imp_type == ImprovementType.BEST_SINCE_FIRST:
            # Compare the (new) best val to the first val
            if self._compute_imp(self.best_val, val) <= 0:
                # No improvement, so use previous best val for comparison
                val = self.best_val
            comp_val = self.first_val
        else:
            raise ValueError(f"Unknown imp_type {self.imp_type}")

        imp = self._compute_imp(comp_val, val)

        return imp

    def get_improvement(self, val: float) -> float:
        imp = self._get_improvement(val)
        if self.transform:
            imp = ImprovementRewardWrapper._pos_log_scale_axis(imp) / (
                np.log10(self.steps + 1) + 1
            )
        return imp

    def step(self, action):
        (
            obs,
            r,
            term,
            trunc,
            info,
        ) = self.env.step(action)
        val = self.get_val(obs, r, info, self.state_eval)
        info[f"raw_{self.state_eval}"] = val
        result_val = self.get_improvement(val)
        self.update_vals(val)
        info[f"best_raw_{self.state_eval}"] = (
            self.best_val if self.best_val is not None else val
        )
        info["improvement"] = result_val
        return obs, result_val, term, trunc, info
