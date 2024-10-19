import gymnasium as gym
import numpy as np
import math
from jaix.env.wrapper import PassthroughWrapper


class AnyFitWrapper(PassthroughWrapper):
    def __init__(self, env: gym.Env, passthrough: bool = True):
        # TODO:could in the future have a config for what exactly to return
        # Could be best target or different scaling etc
        super().__init__(env, passthrough)
        self.steps = 0
        self.best_r = np.inf

    def reset(self, **kwargs):
        self.best_r = np.inf
        self.steps = 0
        self.first_r = None
        return self.env.reset(**kwargs)

    def _pos_log_scale_axis(val):
        if abs(val) <= 1:
            return 0
        ret = math.copysign(1, val) * np.log10(abs(val))
        return ret

    def transform_r(self, r):
        self.best_r = np.min([self.best_r, r])

        imp = self.first_r - self.best_r
        tr = AnyFitWrapper._pos_log_scale_axis(imp) / (np.log10(self.steps) + 1)
        return tr

    def step(self, action):
        (
            obs,
            r,
            term,
            trunc,
            info,
        ) = self.env.step(action)
        if self.steps == 0:
            # Save first reward
            self.first_r = r
            self.best_r = r
        self.steps += 1
        info["raw_r"] = r
        return obs, self.transform_r(r), term, trunc, info
