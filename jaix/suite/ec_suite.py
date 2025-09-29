from jaix.env.singular.ec_env import (
    ECEnvironment,
    ECEnvironmentConfig,
)
from jaix.env.utils.problem.static_problem import StaticProblem
from ttex.config import ConfigurableObjectFactory as COF, Config
from typing import Type, Optional, List
from jaix.suite.suite import Suite, SuiteConfig


class ECSuiteConfig(SuiteConfig):
    def __init__(
        self,
        func_classes: List[Type[StaticProblem]],
        func_configs: List[Config],
        env_config: ECEnvironmentConfig,
        instances: Optional[List[int]] = None,
        agg_instances: Optional[int] = None,
    ):
        self.func_configs = func_configs
        self.func_classes = func_classes
        assert len(func_classes) == len(func_configs)
        functions = list(range(len(func_classes)))
        instances = list(range(15)) if instances is None else instances

        super().__init__(
            env_class=ECEnvironment,
            env_config=env_config,
            functions=functions,
            instances=instances,
            agg_instances=agg_instances,
        )


class ECSuite(Suite):
    config_class = ECSuiteConfig  # type: ignore[assignment]

    def _get_env(self, func, inst):
        func_obj = COF.create(self.func_classes[func], self.func_configs[func], inst)
        return COF.create(ECEnvironment, self.env_config, func_obj, func, inst)
