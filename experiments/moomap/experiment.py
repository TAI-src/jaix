from jaix.env.utils.archive.moomap_archive import MoomapArchive, MoomapArchiveConfig
from jaix.env.wrapper.archive_action_wrapper import (
    ArchiveActionWrapper,
    ArchiveActionWrapperConfig,
)
from jaix.env.wrapper.archive_wrapper import ArchiveWrapper, ArchiveWrapperConfig
from jaix.env.utils.ec_environment import ECEnvironment, ECEnvironmentConfig
from jaix.problem.cobi_problem import CobiProblem, CobiProblemConfig
from ttex.config import ConfigurableObjectFactory as COF, Config
from jaix.env.utils.archive.action_space import ArchiveActionSpace

cobi_config = CobiProblemConfig(
    n_var=3, n_constraints={"Linear": 0, "Quadratic": 0, "Multi": 0}, domain=(-4, 4)
)
ec_env_config = ECEnvironmentConfig(budget_multiplier=1)
archive_config = MoomapArchiveConfig(np_bin=1, coverage_weight=0.5)


class MoomapXConfig(Config):
    def __init__(
        self,
        cobi_config: CobiProblemConfig,
        ec_env_config: ECEnvironmentConfig,
        archive_config: MoomapArchiveConfig,
        action_space_class: type[ArchiveActionSpace],
        action_space_config: Config,
    ):
        super().__init__()
        self.cobi_config = cobi_config
        self.ec_env_config = ec_env_config
        self.archive_config = archive_config
        self.archive_wrapper_config = ArchiveWrapperConfig(
            archive_class=MoomapArchive,
            archive_config=archive_config,
            replace_reward=False,
            passthrough=True,
        )
        self.archive_action_wrapper_config = ArchiveActionWrapperConfig(
            archive_wrapper_config=self.archive_wrapper_config,
            action_space_class=action_space_class,
            action_space_config=action_space_config,
        )


class MoomapX:

    @staticmethod
    def create_env(config: MoomapXConfig):
        func = COF.create(CobiProblem, config.cobi_config, 1)
        env = COF.create(ECEnvironment, config.ec_env_config, func, 0, 1)

        wrapped_env = ArchiveActionWrapper(config.archive_action_wrapper_config, env)

        return wrapped_env

    @staticmethod
    def run(config: MoomapXConfig):
        env = MoomapX.create_env(config)
        obs = env.reset()

        # Prefill the archive with random actions
        for _ in range(10):
            action = env.action_space.sample()
            # TODO: step without wrapper to avoid double translation of action
            obs, reward, terminated, truncated, info = env.step(action)

        # Now measure the success of archive actions
        for _ in range(10):
            action = env.action_space.sample()  # Sample from the archive action space
            # TODO: get current state of the archive
            obs, reward, terminated, truncated, info = env.step(action)
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
