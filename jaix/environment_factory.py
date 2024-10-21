from ttex.config import Config, ConfigurableObjectFactory as COF
from typing import Type, Optional, Union, Dict, Tuple, List
from jaix.suite import Suite, AggType
from jaix.env.composite import CompositeEnvironment
from jaix.env.wrapper import WrappedEnvFactory as WEF, ClosingWrapper, OnlineWrapper
import gymnasium as gym


class CompositeEnvironmentConfig(Config):
    def __init__(
        self,
        agg_type: AggType,
        comp_env_class: Type[CompositeEnvironment],
        comp_env_config: Config,
    ):
        self.agg_type = agg_type
        self.comp_env_class = comp_env_class
        self.comp_env_config = comp_env_config


class EnvironmentConfig(Config):
    default_wrappers = [
        (ClosingWrapper, {}),
        (OnlineWrapper, {"online": True}),
    ]
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
        self.env_wrappers = EnvironmentConfig.default_wrappers + tmp_wrappers

        self.seed = EnvironmentConfig.default_seed if seed is None else seed


class EnvironmentFactory:
    @staticmethod
    def get_envs(env_config: EnvironmentConfig):
        # TODO: potentially add batching here later
        suite = COF.create(env_config.suite_class, env_config.suite_config)
        if env_config.comp_config is None:
            # No composite environments
            for env in suite.get_envs():
                wrapped_env = WEF.wrap(env, env_config.env_wrappers)
                # TODO: reset with seeding here
                yield wrapped_env
                assert wrapped_env.closed
        else:
            comp_config = env_config.comp_config
            for envs in suite.get_agg_envs(
                agg_type=comp_config.agg_type, seed=env_config.seed
            ):
                wrapped_envs = [WEF.wrap(env, env_config.env_wrappers) for env in envs]
                # TODO: The space definitions need to move out of here, this is an assumption
                comp_env = COF.create(
                    comp_config.comp_env_class,
                    comp_config.comp_env_config,
                    wrapped_envs,
                    observation_space=wrapped_envs[0].observation_space,
                    action_space=wrapped_envs[0].action_space,
                )
                wrapped_env = WEF.wrap(comp_env, env_config.env_wrappers)
                # TODO: reset with seeding here
                yield wrapped_env
                assert wrapped_env.closed
                assert all([env.closed for env in wrapped_envs])