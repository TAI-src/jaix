from ttex.config import Config, ConfigurableObjectFactory as COF
from jaix.runner import Runner, Optimiser
from typing import Type
from ttex.log.logging_setup import initiate_logger
from jaix import EnvironmentConfig
from jaix import EnvironmentFactory as EF


class ExperimentConfig(Config):
    def __init__(
        self,
        env_config: EnvironmentConfig,
        runner_class: Type[Runner],
        runner_config: Config,
        opt_class: Type[Optimiser],
        opt_config: Config,
        log_level: int = 30,
    ):
        self.env_config = env_config
        self.runner_class = runner_class
        self.runner_config = runner_config
        self.opt_class = opt_class
        self.opt_config = opt_config
        self.log_level = log_level


class Experiment:
    @staticmethod
    def run(exp_config: ExperimentConfig, *args, **kwargs):
        # Set up logging
        initiate_logger(exp_config.log_level)
        # Log full config here

        runner = COF.create(exp_config.runner_class, exp_config.runner_config)
        for env in EF.get_envs(exp_config.env_config):
            runner.run(
                env, exp_config.opt_class, exp_config.opt_config, *args, **kwargs
            )
            env.close()
