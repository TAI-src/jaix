import gymnasium as gym
from typing import Optional, Dict

class OnlineWrapper(gym.Wrapper):
    def __init__(self, env: gym.Env, online: bool):
        super().__init__(env)
        self.online = online

    def reset(self, options: Optional[Dict], **kwargs):
        reset_opts = {"online": self.online}
        if options is not None:
            reset_opts.update(options)

        return self.env.reset(options=reset_opts, **kwargs)

    def stop(self):
        return self.env.unwrapped.stop()

 




