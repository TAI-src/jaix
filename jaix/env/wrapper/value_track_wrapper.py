from typing import Any

import gymnasium as gym

from jaix.env.wrapper.passthrough_wrapper import PassthroughWrapper


class ValueTrackWrapper(PassthroughWrapper):
    """
    A wrapper that tracks the value of a specified state evaluation key.
    """

    def __init__(
        self,
        env: gym.Env,
        state_eval: list[str] | str = "obs0",
        is_min: dict[str, bool] | bool = True,
        passthrough: bool = True,
    ):
        PassthroughWrapper.__init__(self, env, passthrough)
        if isinstance(state_eval, str):
            self.state_eval = [state_eval]
        else:
            self.state_eval = state_eval.copy()
        if isinstance(is_min, bool):
            self.is_min: dict[str, bool] = {se: is_min for se in self.state_eval}
        else:
            assert is_min.keys() == set(
                self.state_eval
            ), "is_min keys must match state_eval"
            self.is_min = is_min.copy()
        self.best_val: dict[str, float] = {}
        self.first_val: dict[str, float] = {}
        self.last_val: dict[str, float] = {}
        self.steps = 0

    def reset(self, **kwargs):
        self.best_val = {}
        self.steps = 0
        self.first_val = {}
        self.last_val = {}
        return self.env.reset(**kwargs)

    @staticmethod
    def get_vals(
        obs: Any,
        r: float,
        info: dict,
        state_eval: list[str] | str,
    ) -> dict[str, float]:
        if isinstance(state_eval, str):
            state_eval = [state_eval]
        if isinstance(obs, tuple) and len(obs) == 2:
            obs = obs[1]

        values = {}

        for se in state_eval:
            if se == "obs0":
                assert len(obs) == 1
                values[se] = obs[0]
            elif se == "r":
                values[se] = r
            elif se in info:
                values[se] = info[se]
            else:
                raise ValueError(f"Unknown state_eval {se}")

        return values

    def step(self, action):
        (
            obs,
            r,
            term,
            trunc,
            info,
        ) = self.env.step(action)
        vals_dict = self.get_vals(obs, r, info, self.state_eval)
        self.update_vals(vals_dict)
        return obs, r, term, trunc, info

    def update_vals(self, vals_dict: dict[str, float]):
        if self.steps == 0:
            # Save first reward
            self.first_val = vals_dict.copy()
            self.best_val = vals_dict.copy()
            self.last_val = vals_dict.copy()
        assert (
            self.best_val.keys() == vals_dict.keys()
        ), "Keys of best_val and vals_dict must match"
        # Update best and last
        for se in self.state_eval:
            if self.is_min[se]:  # smaller is better
                self.best_val[se] = min(self.best_val[se], vals_dict[se])
            else:  # larger is better
                self.best_val[se] = max(self.best_val[se], vals_dict[se])
        self.last_val = vals_dict.copy()
        self.steps += 1
