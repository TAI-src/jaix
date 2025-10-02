from ttex.config import Config, ConfigurableObjectFactory as COF
from jaix.runner.runner import Runner
from jaix.runner.optimiser import Optimiser
from typing import Type, Optional, Dict
from ttex.log import (
    initiate_logger,
    get_logging_config,
    log_wandb_init,
    teardown_wandb_logger,
)
from jaix.environment_factory import EnvironmentConfig, EnvironmentFactory as EF
import jaix.utils.globals as globals
import logging


class LoggingConfig(Config):
    def __init__(
        self,
        log_level: int = 30,
        logger_name: Optional[str] = None,
        dict_config: Optional[Dict] = None,
    ):
        self.log_level = log_level
        self.logger_name = logger_name if logger_name else globals.LOGGER_NAME
        self.dict_config = (
            dict_config if dict_config else get_logging_config(self.logger_name, False)
        )

    def _setup(self):
        initiate_logger(
            log_level=self.log_level,
            logger_name=self.logger_name,
            disable_existing=False,
            logging_config=self.dict_config,
        )
        return True


class ExperimentConfig(Config):
    def __init__(
        self,
        env_config: EnvironmentConfig,
        runner_class: Type[Runner],
        runner_config: Config,
        opt_class: Type[Optimiser],
        opt_config: Config,
        logging_config: LoggingConfig,
    ):
        self.env_config = env_config
        self.runner_class = runner_class
        self.runner_config = runner_config
        self.opt_class = opt_class
        self.opt_config = opt_config
        self.logging_config = logging_config
        self.run = None

    def setup(self):
        # override to ensure we have a sensible order
        self.logging_config.setup()
        self.env_config.setup()
        self.runner_config.setup()
        self.opt_config.setup()

        # Init wandb if needed
        try:  # TODO: trycatch is tempororary until config._to_dict exists
            config_dict = self.to_dict()
            run = log_wandb_init(
                run_config=config_dict, logger_name=globals.WANDB_LOGGER_NAME
            )
            self.run = run
            if run:
                logging.getLogger(globals.LOGGER_NAME).info(
                    f"Wandb run {run.id} initialized"
                )
            else:
                logging.getLogger(globals.LOGGER_NAME).info("Wandb not initialized")
        except NotImplementedError:
            logging.getLogger(globals.LOGGER_NAME).info(
                "Wandb not installed, skipping wandb logging"
            )
        return True

    def teardown(self):
        self.env_config.teardown()
        self.runner_config.teardown()
        self.opt_config.teardown()
        self.logging_config.teardown()

        teardown_wandb_logger(name=globals.WANDB_LOGGER_NAME)
        return True


class Experiment:
    @staticmethod
    def run(exp_config: ExperimentConfig, *args, **kwargs):
        # Setup experiment
        exp_config.setup()
        logger = logging.getLogger(globals.LOGGER_NAME)

        runner = COF.create(exp_config.runner_class, exp_config.runner_config)
        logger.debug(f"Runner created {runner}")
        for env in EF.get_envs(exp_config.env_config):
            logger.debug(f"Running on env {env}")
            runner.run(
                env, exp_config.opt_class, exp_config.opt_config, *args, **kwargs
            )
            logger.debug(f"Environment {env} done")
            env.close()

        logger.debug("Experiment done")
        exp_config.teardown()
        logger.debug("Experiment torn down")
