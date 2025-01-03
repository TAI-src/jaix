from jaix.env.singular import (
    ECEnvironment,
    ECEnvironmentConfig,
)
from jaix.env.utils.problem import StaticProblem
from ttex.config import ConfigurableObjectFactory as COF, Config
from typing import Type, Optional
from jaix.suite import Suite, AggType


class ECSuiteConfig(Config):
    def __init__(
        self,
        func_class: Type[StaticProblem],
        func_config: Config,
        env_config: ECEnvironmentConfig,
        num_instances: int,
        num_agg_instances: int,
    ):
        self.func_config = func_config
        # TODO: should probably allow multiple functions
        self.env_config = env_config
        self.func_class = func_class
        self.num_instances = num_instances
        self.num_agg_instances = num_agg_instances


class ECSuite(Suite):
    config_class = ECSuiteConfig  # type: ignore[assignment]

    def _get_env(self, inst):
        func = COF.create(self.func_class, self.func_config, inst)
        return COF.create(ECEnvironment, self.env_config, func)
