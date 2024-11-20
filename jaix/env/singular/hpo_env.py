from ttex.config import ConfigurableObject, Config
import gymnasium as gym


class HPOEnvironmentConfig(Config):
    def __init__(self, budget: int, context: str = "D244_F3_C1530_30"):
        self.budget = budget
        self.context = context


class HPOEnvironment(ConfigurableObject, gym.Env):
    config_class = HPOEnvironmentConfig

    def __init__(self, config: HPOEnvironmentConfig, model: str):
        pass
