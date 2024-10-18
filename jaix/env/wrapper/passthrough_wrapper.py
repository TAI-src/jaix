import gymnasium as gym

class PassthroughWrapper(gym.Wrapper):
    def __init__(self, env:gym.Env, passthrough:bool):
        super().__init__(env)
        self.passthrough = passthrough

    def __getattr__(self, name):
        # Pass through anything else that is not overriden
        # This function is only called after it is not found in self
        if self.passthrough: # Passing through
            return getattr(self.env, name)
        else: # Stay consistent with default behaviour
            raise AttributeError(name)
