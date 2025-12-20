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
from uuid import uuid4
from jaix.utils.exp_id import set_exp_id
from jaix.utils.approach_name import set_approach_name
from jaix.runner.ask_tell.at_optimiser import ATOptimiserConfig


class LoggingConfig(Config):
    def __init__(
        self,
        log_level: int = 30,
        logger_name: Optional[str] = None,
        disable_existing: Optional[bool] = False,
        dict_config: Optional[Dict] = None,
    ):
        self.log_level = log_level
        self.logger_name = logger_name if logger_name else globals.LOGGER_NAME
        self.dict_config = (
            dict_config if dict_config else get_logging_config(self.logger_name, False)
        )
        self.disable_existing = disable_existing

    def _setup(self):
        initiate_logger(
            log_level=self.log_level,
            logger_name=self.logger_name,
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
        self.run = None

    def setup(self):
        default_algo_name = f"{self.opt_class.__name__}"
        # TODO: ugly workaround for a good name
        if isinstance(self.opt_config, ATOptimiserConfig):
            default_algo_name = self.opt_config.strategy_class.__name__
        set_approach_name(default_algo_name)

        # override to ensure we have a sensible order
        self.logging_config.setup()
        self.env_config.setup()
        self.runner_config.setup()
        self.opt_config.setup()

        # Init wandb if needed
        try:  # TODO: trycatch is temporary until config._to_dict exists
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
        # override to ensure we have a sensible order
        self.env_config.teardown()
        self.runner_config.teardown()
        self.opt_config.teardown()
        self.logging_config.teardown()

        teardown_wandb_logger(name=globals.WANDB_LOGGER_NAME)
        set_approach_name(None)
        set_exp_id(None)
        return True


class Experiment:
    @staticmethod
    def run(
        exp_config: ExperimentConfig, exp_id: Optional[str] = None, *args, **kwargs
    ):
        # Set up for everything in config, including logging
        exp_config.setup()
        logger = logging.getLogger(globals.LOGGER_NAME)

        # Set experiment ID
        if exp_id is None:
            # If exp_id is not set, first check if we have a wandb run and use that id
            if exp_config.run is not None:
                exp_id = exp_config.run.id
            else:  # Otherwise, generate id
                exp_id = str(uuid4())
        assert exp_id is not None
        set_exp_id(exp_id)

        logger.info(f"Experiment setup with ID {exp_id}")
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
