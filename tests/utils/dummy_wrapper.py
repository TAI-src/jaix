from jaix.env.wrapper import PassthroughWrapper
import gymnasium as gym


class DummyWrapper(PassthroughWrapper):
    def __init__(self, env: gym.Env):
        super().__init__(env, passthrough=True)
        self.stop_dict = {}
        self.env_steps = 0

    def _stop(self):
        return self.stop_dict

    def step(self, action):
        obs, r, term, trunc, info = self.env.step(action)
        self.env_steps += 1
        info["env_step"] = self.env_steps + 1
        return obs, r, term, trunc, info
