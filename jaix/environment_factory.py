from ttex.config import Config, ConfigurableObjectFactory as COF
from typing import Type, Optional, Union, Dict, Tuple, List, Any
from jaix.suite.suite import Suite, AggType
from jaix.env.composite.composite_environment import CompositeEnvironment
from jaix.env.wrapper.wrapped_env_factory import WrappedEnvFactory as WEF
from jaix.env.wrapper.closing_wrapper import ClosingWrapper
import gymnasium as gym
import logging
import jaix.utils.globals as globals

logger = logging.getLogger(globals.LOGGER_NAME)


class CompositeEnvironmentConfig(Config):
    default_wrappers = [
        (ClosingWrapper, {}),
    ]  # type: List[Tuple[Type[gym.Wrapper], Union[Config, Dict]]]

    def __init__(
        self,
        agg_type: AggType,
        comp_env_class: Type[CompositeEnvironment],
        comp_env_config: Config,
        comp_env_wrappers: Optional[
            List[Tuple[Type[gym.Wrapper], Union[Config, Dict[str, Any]]]]
        ] = None,
    ):
        self.agg_type = agg_type
        self.comp_env_class = comp_env_class
        self.comp_env_config = comp_env_config

        # Append default default_wrappers
        tmp_wrappers = [] if comp_env_wrappers is None else comp_env_wrappers
        self.comp_env_wrappers = (
            tmp_wrappers + CompositeEnvironmentConfig.default_wrappers
        )

    def _setup(self):
        success = True
        for _, wrapper_conf in self.comp_env_wrappers:
            if isinstance(wrapper_conf, Config):
                tmp_config: Config = wrapper_conf  # For pyright
                success = tmp_config.setup() and success
        return success

    def _teardown(self):
        success = True
        for _, wrapper_conf in self.comp_env_wrappers:
            if isinstance(wrapper_conf, Config):
                tmp_config: Config = wrapper_conf  # For pyright
                success = tmp_config.teardown() and success
        return success


class EnvironmentConfig(Config):
    default_wrappers = [
        (ClosingWrapper, {}),
    ]  # type: List[Tuple[Type[gym.Wrapper], Union[Config, Dict]]]
    default_seed = 1337

    # TODO: Seeding wrapper
    def __init__(
        self,
        suite_class: Type[Suite],
        suite_config: Config,
        env_wrappers: Optional[
            List[Tuple[Type[gym.Wrapper], Union[Config, Dict]]]
        ] = None,
        comp_config: Optional[CompositeEnvironmentConfig] = None,
        seed: Optional[int] = None,
    ):
        self.suite_class = suite_class
        self.suite_config = suite_config
        self.comp_config = comp_config

        # Append default default_wrappers
        tmp_wrappers = [] if env_wrappers is None else env_wrappers
        self.env_wrappers = tmp_wrappers + EnvironmentConfig.default_wrappers

        self.seed = EnvironmentConfig.default_seed if seed is None else seed

    def _setup(self):
        success = True
        for _, wrapper_conf in self.env_wrappers:
            if isinstance(wrapper_conf, Config):
                tmp_config: Config = wrapper_conf  # For pyright
                success = tmp_config.setup() and success
        return success

    def _teardown(self):
        success = True
        for _, wrapper_conf in self.env_wrappers:
            if isinstance(wrapper_conf, Config):
                tmp_config: Config = wrapper_conf  # For pyright
                success = tmp_config.teardown() and success
        return success


class EnvironmentFactory:
    @staticmethod
    def get_envs(env_config: EnvironmentConfig):
        # TODO: potentially add batching here later
        suite = COF.create(env_config.suite_class, env_config.suite_config)
        logger.debug(f"Suite {suite} created")
        if env_config.comp_config is None:
            # No composite environments
            for env in suite.get_envs():
                logger.debug(f"Environment from suite {env}")
                wrapped_env = WEF.wrap(env, env_config.env_wrappers)
                logger.debug(f"Wrapped env {env}")
                # TODO: reset with seeding here
                yield wrapped_env
                assert wrapped_env.closed
        else:
            comp_config = env_config.comp_config
            for envs in suite.get_agg_envs(
                agg_type=comp_config.agg_type, seed=env_config.seed
            ):
                logger.debug(f"Got {len(envs)} from suite {suite}")
                wrapped_envs = [
                    WEF.wrap(env, comp_config.comp_env_wrappers) for env in envs
                ]
                logger.debug(f"Wrapped envs {wrapped_envs}")
                # Create composite environment
                comp_env = COF.create(
                    comp_config.comp_env_class,
                    comp_config.comp_env_config,
                    wrapped_envs,
                )
                logger.debug(f"Created composite env {comp_env}")
                wrapped_env = WEF.wrap(comp_env, env_config.env_wrappers)
                logger.debug(f"Wrapped composite env {wrapped_env}")
                # TODO: reset with seeding here
                yield wrapped_env
                assert wrapped_env.closed
                assert all([env.closed for env in wrapped_envs])
