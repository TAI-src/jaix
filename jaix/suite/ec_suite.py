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

    def get_envs(self, agg_type: AggType = AggType.NONE, seed: Optional[int] = None):
        if agg_type != AggType.INST and agg_type != AggType.NONE:
            raise NotImplementedError()
        if agg_type == AggType.NONE:
            num_instances = 1
        else:
            num_instances = self.num_instances
        for _ in range(1):
            funcs = [
                COF.create(self.func_class, self.func_config)
                for _ in range(num_instances)
            ]
            envs = [COF.create(ECEnvironment, self.env_config, func) for func in funcs]
            if agg_type == AggType.NONE:
                yield envs[0]
            else:
                yield envs
            assert all([env.closed for env in envs])
