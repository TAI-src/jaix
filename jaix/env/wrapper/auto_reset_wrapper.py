import gymnasium as gym
from ttex.config import ConfigurableObject, Config
from jaix.env.wrapper.passthrough_wrapper import PassthroughWrapper
from typing import Optional


class AutoResetWrapperConfig(Config):
    def __init__(
        self,
        min_steps: int = 1,
        passthrough: bool = True,
        failed_resets_thresh: int = 1,
    ):
        self.min_steps = min_steps
        self.passthrough = passthrough
        self.failed_resets_thresh = failed_resets_thresh
        assert self.failed_resets_thresh > 0


class AutoResetWrapper(PassthroughWrapper, ConfigurableObject):
    """
    A wrapper that automatically resets the environment when done.
    It keeps track of manual and automatic resets, as well as failed resets.
    It also tracks whether the environment was terminated or truncated.
    """

    config_class = AutoResetWrapperConfig

    def __init__(self, config: AutoResetWrapperConfig, env: gym.Env):
        ConfigurableObject.__init__(self, config)
        PassthroughWrapper.__init__(self, env, self.passthrough)
        self.man_resets = 0
        self.auto_resets = 0
        self.failed_resets = 0
        self.steps = 0
        self.prev_r = None
        self.trunc = False
        self.term = False

    def reset(
        self,
        *,
        seed: Optional[int] = None,
        options: Optional[dict] = None,
    ):
        obs, info = self.env.reset(seed=seed, options=options)
        if options is None or "auto" not in options or not options["auto"]:
            # Manual reset
            self.man_resets += 1
            if self.prev_r is not None:
                info["final_r"] = self.prev_r
        else:
            # Auto reset
            self.auto_resets += 1
            if self.steps >= self.min_steps:
                # Only update final r if it is not a failed reset
                info["final_r"] = self.prev_r
        self.steps = 0
        self.prev_r = None
        if options is None or "online" not in options or not options["online"]:
            # RL reset, full reset of term and trunc
            self.trunc = False
            self.term = False
        else:
            # Online reset, so only termination (success), not truncation due to time constraints
            self.term = False
        return obs, info

    def step(self, action):
        (
            obs,
            r,
            term,
            trunc,
            info,
        ) = self.env.step(action)
        if term:
            # from https://gymnasium.farama.org/_modules/gymnasium/wrappers/autoreset/
            _, reset_info = self.reset(options={"online": True, "auto": True})
            assert (
                "final_observation" not in reset_info
            ), 'info dict cannot contain key "final_observation" '
            assert (
                "final_info" not in reset_info
            ), 'info dict cannot contain key "final_info" '
            info_updates = {"final_observation": obs, "final_info": info}
            if "final_r" in reset_info:
                info_updates["final_r"] = reset_info["final_r"]

            (
                obs,
                r,
                term,
                trunc,
                info,
            ) = self.env.step(self.action_space.sample())

            info.update(info_updates)

            if "final_r" not in info:
                # This means we reset previously and on the first step
                # we are done. That is a fail
                self.failed_resets += 1
        elif trunc:
            info["final_observation"] = obs
            info["final_info"] = info
            info["final_r"] = r

        self.steps += 1
        self.prev_r = r
        if term:
            self.term = True
        if trunc:
            self.trunc = True
        return obs, r, term, trunc, info

    def _stop(self):
        stop_dict = {}
        if self.failed_resets >= self.failed_resets_thresh:
            stop_dict["failed_resets"] = self.failed_resets
        if self.term:
            stop_dict["term"] = True
        if self.trunc:
            stop_dict["trunc"] = True
        return stop_dict
