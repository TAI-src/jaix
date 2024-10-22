from jaix.env.wrapper import PassthroughWrapper
import gymnasium as gym


class DummyWrapper(PassthroughWrapper):
    def __init__(self, env: gym.Env):
        super().__init__(env, passthrough=True)
        self.stop_dict = {}

    def _stop(self):
        return self.stop_dict
