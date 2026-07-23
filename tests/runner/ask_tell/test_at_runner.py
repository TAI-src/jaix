import pytest
from jaix.env.wrapper.wrapped_env_factory import WrappedEnvFactory as WEF
from jaix.environment_factory import EnvironmentConfig as EC
from jaix.runner.ask_tell.ask_tell_runner import ATRunner, ATRunnerConfig
from jaix.runner.ask_tell.at_optimiser import ATOptimiser, ATOptimiserConfig
from jaix.runner.ask_tell.strategy.cma import CMA, CMAConfig
from jaix.runner.ask_tell.strategy.random import (
    RandomATStrat,
    RandomATStratConfig,
)
from ttex.config import ConfigurableObjectFactory as COF

from . import DummyEnv


def get_optimiser(opts: str = None):
    if opts == "CMA":
        config = ATOptimiserConfig(
            strategy_class=CMA,
            strategy_config=CMAConfig(sigma0=5),
            init_pop_size=1,
            stop_after=3,
        )
    else:
        config = ATOptimiserConfig(
            strategy_class=RandomATStrat,
            strategy_config=RandomATStratConfig(ask_size=5),
            init_pop_size=1,
            stop_after=3,
        )
    return config


@pytest.mark.parametrize("opts", ["CMA", "Random"])
def test_run(opts):
    wrappers = EC.default_wrappers

    env = WEF.wrap(DummyEnv(), wrappers)
    opt_config = get_optimiser(opts)
    runner = COF.create(ATRunner, ATRunnerConfig(max_evals=4, disp_interval=50))
    runner.run(env, ATOptimiser, opt_config)
