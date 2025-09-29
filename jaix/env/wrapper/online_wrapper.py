import gymnasium as gym
from typing import Optional, Dict
from jaix.env.wrapper.passthrough_wrapper import PassthroughWrapper


class OnlineWrapper(PassthroughWrapper):
    def __init__(self, env: gym.Env, online: bool, passthrough: bool = True):
        super().__init__(env, passthrough)
        self.online = online

    def reset(self, *, seed: Optional[int] = None, options: Optional[Dict] = None):
        reset_opts = {"online": self.online}
        if options is not None:
            reset_opts.update(options)

        return self.env.reset(seed=seed, options=reset_opts)
