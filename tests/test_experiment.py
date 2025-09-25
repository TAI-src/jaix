from .test_environment_factory import comp_config, ec_config, env_config
from .runner.ask_tell.test_at_runner import get_optimiser
import pytest
from jaix.runner.ask_tell import ATRunnerConfig, ATRunner, ATOptimiser
from jaix import ExperimentConfig, Experiment, LoggingConfig


def exp_config(ec_config, comp_config, comp: bool, opts: str = None):
    if comp:
        env_conf = env_config(ec_config, comp_config=comp_config)
    else:
        env_conf = env_config(ec_config)
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


@pytest.mark.parametrize("comp", [False, True])
def test_experiment(ec_config, comp_config, comp):
    config = exp_config(ec_config, comp_config, comp=comp, opts="Random")
    Experiment.run(config)
