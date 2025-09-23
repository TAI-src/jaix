from ttex.config import Config, ConfigurableObjectFactory as COF
from jaix.runner import Runner, Optimiser
from typing import Type, Optional, Dict
from ttex.log import initiate_logger, get_logging_config
from jaix import EnvironmentConfig, EnvironmentFactory as EF, LOGGER_NAME
import logging


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
    def run(exp_config: ExperimentConfig, *args, **kwargs):
        # Set up for everything in config, including logging
        exp_config.setup()
        print("HELLOOO________________________----")
        print(exp_config.env_config.env_wrappers[1][1].__dict__)
        logger = logging.getLogger(LOGGER_NAME)

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
        return True
