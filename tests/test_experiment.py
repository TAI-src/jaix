from .test_environment_factory import comp_config, ec_config, env_config
from .runner.ask_tell.test_at_runner import get_optimiser
import pytest
from jaix.runner.ask_tell.ask_tell_runner import ATRunnerConfig, ATRunner
from jaix.runner.ask_tell.at_optimiser import ATOptimiser
from jaix.experiment import ExperimentConfig, Experiment
from jaix.utils.logging_config import LoggingConfig
import shutil
from .utils.dummy_wrapper import DummyWrapper, DummyWrapperConfig
from jaix.environment_factory import EnvironmentFactory as EF
from jaix.env.wrapper.closing_wrapper import ClosingWrapper
from jaix.utils.experiment_context import ExperimentContext


def exp_config(
    ec_config, comp_config, comp: bool, opts: str = "Random"
) -> ExperimentConfig:
    wrappers = [
        (DummyWrapper, DummyWrapperConfig(passthrough=True)),
    ]
    if comp:
        env_conf = env_config(ec_config, comp_config=comp_config, wrappers=wrappers)
    else:
        env_conf = env_config(ec_config, wrappers=wrappers)
    opt_config = get_optimiser(opts)
    runner_config = ATRunnerConfig(max_evals=4, disp_interval=50)
    config = ExperimentConfig(
        env_config=env_conf,
        runner_class=ATRunner,
        runner_config=runner_config,
        opt_class=ATOptimiser,
        opt_config=opt_config,
        logging_config=LoggingConfig(log_level=10),
    )
    return config


def test_experiment_config(ec_config, comp_config):
    ctx = ExperimentContext()
    config = exp_config(ec_config, comp_config, comp=False, opts="Random")
    assert ctx.get("approach_name") is None

    assert config.setup(ctx)
    assert config.env_config.env_wrappers[0][1]._stp is True

    default_algo_name = ctx.get("approach_name")
    assert default_algo_name is not None

    assert config.teardown(ctx)
    assert config.env_config.env_wrappers[0][1].trdwn is True


@pytest.mark.parametrize("comp", [False, True])
def test_experiment(ec_config, comp_config, comp):
    config = exp_config(ec_config, comp_config, comp=comp, opts="Random")
    exp_id = Experiment.run(config)
    assert exp_id is not None
