from jaix.env.singular import (
    ECEnvironment,
    ECEnvironmentConfig,
)
from jaix.env.utils.problem import StaticProblem
from ttex.config import ConfigurableObject, ConfigurableObjectFactory as COF, Config
from typing import Type, Optional
from jaix.suite import Suite, AggType


class ECSuiteConfig(Config):
    def __init__(
        self,
        func_class: Type[StaticProblem],
        func_config: Config,
        env_config: ECEnvironmentConfig,
        num_instances: int = 1,
    ):
        self.func_config = func_config
        self.env_config = env_config
        self.func_class = func_class
        self.num_instances = num_instances


class ECSuite(ConfigurableObject, Suite):
    config_class = ECSuiteConfig

    def get_envs(self):
        for _ in range(1):
            func = COF.create(self.func_class, self.func_config)
            env = COF.create(ECEnvironment, self.env_config, func)
            yield env

    def get_agg_envs(self, agg_type: AggType, seed: Optional[int] = None):
        for _ in range(1):
            funcs = [
                COF.create(self.func_class, self.func_config)
                for _ in range(self.num_instances)
            ]
            envs = [COF.create(ECEnvironment, self.env_config, func) for func in funcs]
            yield envs
