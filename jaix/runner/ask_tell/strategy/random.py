import numpy as np
from jaix.runner.ask_tell.at_strategy import ATStrategy
from ttex.config import Config, ConfigurableObject


class RandomATStratConfig(Config):
    def __init__(self, ask_size: int):
        self.ask_size = ask_size


class RandomATStrat(ConfigurableObject, ATStrategy):
    config_class = RandomATStratConfig

    def __init__(self, config: RandomATStratConfig, xstart, *args, **kwargs):
        ConfigurableObject.__init__(self, config)
        ATStrategy.__init__(self, xstart)

    def initialize(self):
        self.xcurrent = [np.array(xi) for xi in self.xstart]

    @property
    def name(self):
        return "Random Search"

    def ask(self, env, **kwargs):
        self.xcurrent = [env.action_space.sample() for _ in range(self.ask_size)]
        return self.xcurrent

    def tell(self, **kwargs):
        pass

    def stop(self):
        return {}
