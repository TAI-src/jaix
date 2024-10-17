import gymnasium as gym
from ttex.config import ConfigurableObject, Config

class AutoResetWrapperConfig(Config):
    def __init__(self, min_steps: int = 1):
        self.min_steps = min_steps

class AutoResetWrapper(gym.Wrapper, ConfigurableObject):
    config_class = AutoResetWrapperConfig

    def __init__(self, config: AutoResetWrapperConfig, env: gym.Env):
        ConfigurableObject.__init__(self, config)
        gym.Wrapper.__init__(self, env)
        self.man_resets = 0
        self.auto_resets = 0
        self.failed_resets = 0
        self.steps = 0
        self.prev_r = None

    def reset(self, **kwargs):
        obs, info = self.env.reset(**kwargs)
        if self.steps >= self.min_steps:
            # Only update final r if it is not a failed reset
            info["final_r"] = self.prev_r
        self.man_resets += 1
        self.steps = 0
        self.prev_r = None
        return obs, info

    def step(self, action):
        (
            obs,
            r,
            term,
            trunc,
            info,
        ) = self.env.step(action)
        if term or trunc:
            # from https://gymnasium.farama.org/_modules/gymnasium/wrappers/autoreset/
            new_obs, new_info = self.env.reset(options={"online": True})
            assert (
                "final_observation" not in new_info
            ), 'info dict cannot contain key "final_observation" '
            assert (
                "final_info" not in new_info
            ), 'info dict cannot contain key "final_info" '
            assert "final_r" not in new_info, 'info dict cannot contain key "final_r" '

            new_info["final_observation"] = obs
            new_info["final_info"] = info

            obs = new_obs
            info = new_info

            if self.steps == 0:
                # This means we reset previously and on the first step
                # we are done. That is a fail
                self.failed_resets += 1
            elif self.steps >= self.min_steps:
                # Only update final r if it is not a failed reset
                info["final_r"] = r
            self.steps = 0
            self.auto_resets += 1
        else:
            self.steps += 1
        self.prev_r = r
        return obs, r, term, trunc, info

    def stop(self, failed_resets_thresh: int = 1):
        # Consider finally stopped if more failed resets than desired
        return self.failed_resets >= failed_resets_thresh
