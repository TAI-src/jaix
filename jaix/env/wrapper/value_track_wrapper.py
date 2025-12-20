from jaix.env.wrapper.passthrough_wrapper import PassthroughWrapper
from typing import Optional, Any, Dict
import gymnasium as gym


class ValueTrackWrapper(PassthroughWrapper):
    """
    A wrapper that tracks the value of a specified state evaluation key.
    """

    def __init__(
        self,
        env: gym.Env,
        state_eval: str = "obs0",
        is_min: bool = True,
        passthrough: bool = True,
    ):
        PassthroughWrapper.__init__(self, env, passthrough)
        self.state_eval = state_eval
        self.is_min = is_min
        self.best_val: Optional[float] = None
        self.steps = 0
        self.first_val: Optional[float] = None
        self.last_val: Optional[float] = None

    def reset(self, **kwargs):
        self.best_val = None
        self.steps = 0
        self.first_val = None
        self.last_val = None
        return self.env.reset(**kwargs)

    @staticmethod
    def get_val(obs: Any, r: float, info: Dict, state_eval: str) -> float:
        if state_eval == "obs0":
            if isinstance(obs, tuple) and len(obs) == 2:
                # This returns a tuple with the current env number and the observation
                # so we need to extract the observation
                # (env_num, ob) = obs
                obs = obs[1]
            assert len(obs) == 1
            val = obs[0]
        elif state_eval == "r":
            val = r
        elif state_eval in info:
            val = info[state_eval]
        else:
            raise ValueError(f"Unknown state_eval {state_eval}")
        return val

    def step(self, action):
        (
            obs,
            r,
            term,
            trunc,
            info,
        ) = self.env.step(action)
        val = self.get_val(obs, r, info, self.state_eval)
        self.update_vals(val)
        return obs, r, term, trunc, info

    def update_vals(self, val):
        if self.steps == 0:
            # Save first reward
            self.first_val = val
            self.best_val = val
            self.last_val = val
        assert self.best_val is not None
        # Update best and last
        self.best_val = (
            min(self.best_val, val) if self.is_min else max(self.best_val, val)
        )
        self.last_val = val
        self.steps += 1
