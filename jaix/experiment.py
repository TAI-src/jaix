from ttex.config import Config, ConfigurableObjectFactory as COF
from jaix.runner.runner import Runner
from jaix.runner.optimiser import Optimiser
from typing import Type, Optional, Dict
from ttex.log import initiate_logger, get_logging_config
from jaix.environment_factory import EnvironmentConfig, EnvironmentFactory as EF
from jaix.utils.globals import LOGGER_NAME
import logging
from uuid import uuid4
from jaix.utils import set_exp_id


class LoggingConfig(Config):
    def __init__(
        self,
        log_level: int = 30,
        logger_name: Optional[str] = None,
        disable_existing: Optional[bool] = False,
        dict_config: Optional[Dict] = None,
    ):
        self.log_level = log_level
        self.disable_existing = disable_existing
        self.logger_name = logger_name if logger_name else LOGGER_NAME
        self.dict_config = (
            dict_config
            if dict_config
            else get_logging_config(self.logger_name, self.disable_existing)
        )

    def _setup(self):
        initiate_logger(
            log_level=self.log_level,
            logger_name=LOGGER_NAME,
            disable_existing=self.disable_existing,
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


class Experiment:
    @staticmethod
    def run(
        exp_config: ExperimentConfig, exp_id: Optional[str] = None, *args, **kwargs
    ):
        # Set experiment ID
        exp_id = exp_id if exp_id is not None else str(uuid4())
        set_exp_id(exp_id)

        # Set up for everything in config, including logging
        exp_config.setup()
        logger = logging.getLogger(LOGGER_NAME)

        logger.error(f"Experiment setup with ID {exp_id}")
        runner = COF.create(exp_config.runner_class, exp_config.runner_config)
        logger.debug(f"Runner created {runner}")
        for env in EF.get_envs(exp_config.env_config):
            logger.debug(f"Running on env {env}")
            if env.stop():
                logger.warning(f"Environment {env} already stopped, skipping")
                env.close()
                continue
            runner.run(
                env, exp_config.opt_class, exp_config.opt_config, *args, **kwargs
            )
            logger.debug(f"Environment {env} done")
            env.close()
        logger.debug(f"Experiment {exp_id} done")

        exp_config.teardown()
        logger.debug(f"Experiment {exp_id} torn down")
        return exp_id
