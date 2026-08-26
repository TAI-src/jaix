from re import A
from jaix.env.utils.archive.moomap_archive import MoomapArchive, MoomapArchiveConfig
from jaix.env.wrapper.archive_action_wrapper import (
    ArchiveActionWrapper,
    ArchiveActionWrapperConfig,
)
from jaix.suite.ec_suite import ECSuite, ECSuiteConfig
from jaix.environment_factory import EnvironmentConfig
from jaix.env.wrapper.archive_wrapper import ArchiveWrapper, ArchiveWrapperConfig
from jaix.env.singular.ec_env import ECEnvironmentConfig
from ttex.config import ConfigurableObjectFactory as COF, Config


from jaix.environment_factory import EnvironmentFactory as EF
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from jaix.env.utils.archive.action_space import (
    UniformCrossoverActionSpace,
    UniformCrossoverActionSpaceConfig,
)
from jaix.env.wrapper.wrapped_env_factory import WrappedEnvFactory as WEF

from ttex.config.config import ConfigFactory as CF
import argparse
import json


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
        self.num_parents = 2  # Number of parents for the uniform crossover action space

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
                crossover_attribute="info.env_action", num_parents=self.num_parents
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
                func_configs=[REProblemConfig()],
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
    def prefill_archive(env, config: MoomapXConfig):
        archive_wrapper_config = config.generate_archive_wrapper_config()
        wrappers = [(ArchiveWrapper, archive_wrapper_config)]
        wa_env = WEF.wrap(env, wrappers)
        for _ in range(config.num_samples):
            action = wa_env.action_space.sample()
            obs, reward, terminated, truncated, info = wa_env.step(action)
        MoomapX.process_archive(wa_env.archive)
        return wa_env.archive

    @staticmethod
    def process_archive(archive: MoomapArchive):
        for entry in archive.get_all():
            # Add original action info to assure same format for entries added from archive action
            info = getattr(entry, "info", {})
            info["env_action"] = getattr(entry, "action", None)
            setattr(entry, "info", info)

    @staticmethod
    def count_successes(archive: MoomapArchive):
        archive_stats = archive.get_archive_stats()
        return (
            archive_stats[f"add_counter_{archive.n_bins}"]
            + archive_stats[f"replace_counter_{archive.n_bins}"]
        )

    @staticmethod
    def run(config: MoomapXConfig):
        env_config = config.generate_env_config()
        for env in EF.get_envs(env_config):
            archive = MoomapX.prefill_archive(env, config)

            archive_action_wrapper_config = (
                config.generate_archive_action_wrapper_config()
            )
            wrappers = [(ArchiveActionWrapper, archive_action_wrapper_config)]
            waa_env = WEF.wrap(env, wrappers)

            waa_env.set_archive(
                archive
            )  # Set the prefilled archive in the wrapped environment

            prev_n_success = MoomapX.count_successes(archive)
            # Now measure the success of archive actions
            for _ in range(config.num_trials):
                action = (
                    waa_env.action_space.sample()
                )  # Sample from the archive action space
                archive_stats = archive.get_archive_stats()
                # TODO: get current state of the archive
                try:
                    obs, reward, terminated, truncated, info = waa_env.step(action)
                except AssertionError as e:
                    assert str(e) == "Not enough parents found in the archive"
                    continue
                new_n_success = MoomapX.count_successes(archive)
                if new_n_success > prev_n_success:
                    print(f"Action {action} was successful!")
                else:
                    print(f"Action {action} was NOT successful. ")
                prev_n_success = new_n_success

            env.close()

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
        #


def parse_args():
    parser = argparse.ArgumentParser(description="Run MoomapX experiment")
    parser.add_argument(
        "--config_file", type=str, help="Path to the configuration file"
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    with open(args.config_file, "r") as f:
        run_config = json.load(f)
    config = CF.from_dict(run_config)
    MoomapX.run(config)
