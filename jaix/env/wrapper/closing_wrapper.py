import gymnasium as gym

class ClosingWrapper(gym.Wrapper):
    
    def __init__(self, env:gym.Env):
        super().__init__(env)
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

