from ttex.config import Config, ConfigurableObjectFactory as COF
from jaix.runner.runner import Runner
from jaix.runner.optimiser import Optimiser
from typing import Type, Optional, Dict
from jaix.environment_factory import EnvironmentConfig, EnvironmentFactory as EF
import logging
from jaix.runner.ask_tell.at_optimiser import ATOptimiserConfig
from jaix.utils.experiment_context import ExperimentContext
from jaix.utils.logging_config import LoggingConfig
from jaix.utils.wandb_session import WandbSession


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
        Config.__init__(self)
        self.env_config = env_config
        self.runner_class = runner_class
        self.runner_config = runner_config
        self.opt_class = opt_class
        self.opt_config = opt_config
        self.logging_config = logging_config
        self.run = None

    def gen_approach_name(self) -> str:
        default_algo_name = f"{self.opt_class.__name__}"
        # TODO: ugly workaround for a good name
        if isinstance(self.opt_config, ATOptimiserConfig):
            default_algo_name = self.opt_config.strategy_class.__name__
        return default_algo_name

    def setup(self, ctx: ExperimentContext):
        ctx.set("approach_name", self.gen_approach_name())

        # override to ensure we have a sensible order
        self.logging_config.setup(ctx)
        self.env_config.setup(ctx)
        self.runner_config.setup(ctx)
        self.opt_config.setup(ctx)

        return True

    def teardown(self, ctx: ExperimentContext):
        # override to ensure we have a sensible order
        self.env_config.teardown(ctx)
        self.runner_config.teardown(ctx)
        self.opt_config.teardown(ctx)
        self.logging_config.teardown(ctx)

        return True


class Experiment:
    @staticmethod
    def run(
        exp_config: ExperimentConfig, exp_id: Optional[str] = None, *args, **kwargs
    ):
        ctx = ExperimentContext(exp_id)

        try:
            # Set up for everything in config, including logging
            exp_config.setup(ctx)
            # Start wandb session if needed
            if ctx.get("wandb_logger_name") is not None:
                wandb_session = WandbSession(exp_config)
                wandb_session.start(ctx)
            ctx.freeze()
            exp_config.set_context(ctx)
            Experiment._run(exp_config, ctx, *args, **kwargs)
        finally:
            exp_config.teardown(ctx)
            # End wandb session if needed
            if ctx.get("wandb_logger_name") is not None:
                wandb_session.end(ctx)
        return ctx.get("exp_id")

    @staticmethod
    def _run(exp_config: ExperimentConfig, ctx: ExperimentContext, *args, **kwargs):
        logger = logging.getLogger(ctx.get("logger_name"))
        logger.info(f"Experiment setup with ID {ctx.get('exp_id')}")
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
        logger.debug(f"Experiment {ctx.get("exp_id")} done")
