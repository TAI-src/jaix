from ttex.config import Config
from ttex.config import ConfigurableObjectFactory as COF

from jaix.env.singular.ec_env import (
    ECEnvironment,
    ECEnvironmentConfig,
)
from jaix.env.utils.problem.static_problem import StaticProblem
from jaix.suite.suite import Suite, SuiteConfig


class ECSuiteConfig(SuiteConfig):
    def __init__(
        self,
        func_classes: list[type[StaticProblem]],
        func_configs: list[Config],
        env_config: ECEnvironmentConfig,
        instances: list[int] | None = None,
        agg_instances: int | None = None,
    ):
        self.func_configs = func_configs
        self.func_classes = func_classes
        assert len(func_classes) == len(func_configs)
        functions = list(range(len(func_classes)))
        instances = list(range(15)) if instances is None else instances

        _agg_instances: list[tuple[int, ...]] | None | int
        if agg_instances is not None and agg_instances == 0:
            _agg_instances = [
                (0,)
            ]  # FIXME: This is a hack to avoid double-initialisation which causes issues when there is only one instance per function. This should be fixed properly in the future.
        else:
            _agg_instances = agg_instances
        super().__init__(
            env_class=ECEnvironment,
            env_config=env_config,
            functions=functions,
            instances=instances,
            agg_instances=_agg_instances,
        )


class ECSuite(Suite):
    config_class = ECSuiteConfig  # type: ignore[assignment]

    def _get_env(self, func, inst):
        func_obj = COF.create(self.func_classes[func], self.func_configs[func], inst)
        return COF.create(ECEnvironment, self.env_config, func_obj, func, inst)
