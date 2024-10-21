import gymnasium as gym
from ttex.config import ConfigurableObject, Config
from jaix.env.wrapper import PassthroughWrapper


class MaxEvalWrapperConfig(Config):
    def __init__(self, max_evals: int, passthrough: bool = True):
        self.max_evals = max_evals
        self.passthrough = passthrough


class AutoResetWrapper(PassthroughWrapper, ConfigurableObject):
    config_class = MaxEvalWrapperConfig

    def __init__(self, config:MaxEvalWrapperConfig, env: gym.Env):
        ConfigurableObject.__init__(self, config)
        PassthroughWrapper.__init__(self, env, self.passthrough)
        self.evals = 0

    def reset(self, **kwargs):
    # TODO: online vs offline reset
    #

    def step(self, action):
        self.evals +=1
        return self.env.step(action)

    def stop(self):
        # TODO: not always present
        term_conds = {}
        if self.countiter >= self.stop_after:
            term_conds["countiter"] = self.countiter
        term_conds.update(self.strategy.stop())
        return term_conds  

