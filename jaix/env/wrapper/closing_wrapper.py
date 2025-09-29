import gymnasium as gym
from jaix.env.wrapper.passthrough_wrapper import PassthroughWrapper


class ClosingWrapper(PassthroughWrapper):
    def __init__(self, env: gym.Env, passthrough: bool = True):
        super().__init__(env, passthrough)
        self.closed = False

    def reset(self, **kwargs):
        if self.closed:
            raise ValueError(f"Environment {self.unwrapped} already closed")
        else:
            return self.env.reset(**kwargs)

    def step(self, *args, **kwargs):
        if self.closed:
            raise ValueError(f"Environment {self.unwrapped} already closed")
        else:
            return self.env.step(*args, **kwargs)

    def close(self, **kwargs):
        res = self.env.close(**kwargs)
        self.closed = True
        return res
