from jaix.env.utils.archive.moomap_archive import MoomapArchive, MoomapArchiveConfig
from jaix.env.wrapper.archive_action_wrapper import (
    ArchiveActionWrapper,
    ArchiveActionWrapperConfig,
)
from jaix.suite.ec_suite import ECSuite, ECSuiteConfig
from jaix.environment_factory import EnvironmentConfig
from jaix.env.wrapper.archive_wrapper import ArchiveWrapper, ArchiveWrapperConfig
from jaix.env.utils.ec_environment import ECEnvironmentConfig
from ttex.config import ConfigurableObjectFactory as COF, Config


from jaix.environment_factory import EnvironmentFactory as EF
from jaix.problem.reproblem import REProblem, ReProblemConfig
from jaix.env.utils.archive.action_space import (
    UniformCrossoverActionSpace,
    UniformCrossoverActionSpaceConfig,
)
from jaix.env.wrapper.wrapped_env_factory import WrappedEnvFactory as WEF


class MoomapXConfig(Config):
    def __init__(
        self,
        moomap_config: MoomapArchiveConfig,  # Configuration for the Moomap archive
        env_budget_multiplier: int,  # Budget multiplier for the EC environment
        num_samples: int,  # Number of random samples to prefill the archive with
        num_trials: int,  # Number of trials to run for measuring the success of archive actions
        mode: str = "reproblem",  # Mode for the environment, either "cobi" or "reproblem"
        seed: int | None = None,  # Seed for reproducibility
    ):
        super().__init__()
        self.moomap_config = moomap_config
        self.env_budget_multiplier = env_budget_multiplier
        self.num_samples = num_samples
        self.num_trials = num_trials
        self.seed = seed
        self.mode = mode

    def generate_archive_wrapper_config(self) -> ArchiveWrapperConfig:
        return ArchiveWrapperConfig(
            archive_class=MoomapArchive,
            archive_config=self.moomap_config,
            replace_reward=False,
            passthrough=True,
        )

    def generate_archive_action_wrapper_config(self) -> ArchiveActionWrapperConfig:
        return ArchiveActionWrapperConfig(
            archive_wrapper_config=self.generate_archive_wrapper_config(),
            action_space_class=UniformCrossoverActionSpace,
            action_space_config=UniformCrossoverActionSpaceConfig(
                crossover_attribute="info.original_action", num_parents=2
            ),
        )

    def generate_env_config(self) -> EnvironmentConfig:

        ec_env_config = ECEnvironmentConfig(
            budget_multiplier=self.env_budget_multiplier
        )
        if self.mode == "cobi":
            raise NotImplementedError("Cobi mode is not implemented yet.")
        elif self.mode == "reproblem":
            ec_suite_config = ECSuiteConfig(
                func_classes=[REProblem],
                func_configs=ReProblemConfig(),
                env_config=ec_env_config,
                instances=list(range(23)),  # 23 instances of ReProblem
                agg_instances=1,
            )

        env_config = EnvironmentConfig(
            suite_class=ECSuite,
            suite_config=ec_suite_config,
            env_wrappers=None,
            comp_config=None,
            seed=self.seed,
        )
        return env_config


class MoomapX:

    @staticmethod
    def run(config: MoomapXConfig):
        env_config = config.generate_env_config()
        for env in EF.get_envs(env_config):
            archive_wrapper_config = config.generate_archive_wrapper_config()
            wa_env = WEF.wrap(env, archive_wrapper_config)
            obs = wa_env.reset()
            for _ in range(config.num_samples):
                action = wa_env.action_space.sample()
                obs, reward, terminated, truncated, info = wa_env.step(action)

            archive_action_wrapper_config = (
                config.generate_archive_action_wrapper_config()
            )
            waa_env = WEF.wrap(env, archive_action_wrapper_config)
            # Transfer archive from wa_env to waa_env
            waa_env.archive = wa_env.archive

            # Now measure the success of archive actions
            for _ in range(config.num_trials):
                action = (
                    waa_env.action_space.sample()
                )  # Sample from the archive action space
                # TODO: get current state of the archive
                obs, reward, terminated, truncated, info = waa_env.step(action)
                # TODO: evaluate the success of the action based on if a new sample was added
                # TODO: also determine the distance in objective and search space between the parents
                #
                # For now, just print the info

        # TODO: turn into a dataframe
        # TODO: figure out which cobi problems are interesting
        # Might still need to create a suite so that I don't have to create a new config for each new one
        #
        #
        #
