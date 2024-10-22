from jaix.runner.ask_tell.strategy import (
    CMAConfig,
    CMA,
    RandomATStrat,
    RandomATStratConfig,
)
from . import DummyEnv
from jaix.runner.ask_tell import ATRunnerConfig, ATRunner
import pytest
from ttex.config import ConfigurableObjectFactory as COF
from jaix.env.wrapper import (
    WrappedEnvFactory as WEF,
    MaxEvalWrapper,
    MaxEvalWrapperConfig,
)
from jaix import EnvironmentConfig as EC
from jaix.runner.ask_tell import ATOptimiserConfig, ATOptimiser


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
