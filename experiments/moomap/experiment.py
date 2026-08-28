import argparse
import json
import os
import uuid
from typing import Any

import numpy as np
import pandas as pd
from jaix.env.singular.ec_env import ECEnvironmentConfig
from jaix.env.utils.archive.action_space import (
    UniformCrossoverActionSpace,
    UniformCrossoverActionSpaceConfig,
)
from jaix.env.utils.archive.moomap_archive import (
    MoomapArchive,
    MoomapArchiveConfig,
    MoomapArchiveEntry,
)
from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from jaix.env.wrapper.archive_action_wrapper import (
    ArchiveActionWrapper,
    ArchiveActionWrapperConfig,
)
from jaix.env.wrapper.archive_wrapper import ArchiveWrapper, ArchiveWrapperConfig
from jaix.env.wrapper.wrapped_env_factory import WrappedEnvFactory as WEF
from jaix.environment_factory import EnvironmentConfig
from jaix.environment_factory import EnvironmentFactory as EF
from jaix.suite.ec_suite import ECSuite, ECSuiteConfig
from ttex.config import Config
from ttex.config.config import ConfigFactory as CF

from cobi_config_generator import get_configs as get_cobi_configs


class MoomapXConfig(Config):
    def __init__(
        self,
        moomap_config: MoomapArchiveConfig,  # Configuration for the Moomap archive
        num_samples: int,  # Number of random samples to prefill the archive with
        num_trials: int,  # Number of trials to run for measuring the success of archive actions
        mode: str = "reproblem",  # Mode for the environment, either "cobi" or "reproblem"
        seed: int | None = None,  # Seed for reproducibility
    ):
        super().__init__()
        self.moomap_config = moomap_config
        self.num_samples = num_samples
        self.num_trials = num_trials
        self.seed = seed
        self.mode = mode
        self.num_parents = 2  # Number of parents for the uniform crossover action space

    def generate_archive_wrapper_config(self) -> ArchiveWrapperConfig:
        return ArchiveWrapperConfig(
            archive_class=MoomapArchive,
            archive_config=self.moomap_config,
            replace_reward=True,
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
            budget_multiplier=1  # This is overriden anyway since we do not stop based on the enviornment
        )
        if self.mode == "cobi":
            cobi_configs = get_cobi_configs()
            ec_suite_config = ECSuiteConfig(
                func_classes=[CobiProblem] * len(cobi_configs),
                func_configs=cobi_configs,
                env_config=ec_env_config,
                instances=[1],  # One instance per Cobi config
                agg_instances=0,
            )
        elif self.mode == "reproblem":
            ec_suite_config = ECSuiteConfig(
                func_classes=[REProblem],
                func_configs=[REProblemConfig()],
                env_config=ec_env_config,
                instances=list(range(23)),  # 23 instances of ReProblem
                agg_instances=0,
            )
        else:
            raise ValueError(f"Unknown mode: {self.mode}")

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
            _obs, _reward, _terminated, _truncated, _info = wa_env.step(action)
        MoomapX.process_archive(wa_env.archive)
        return wa_env.archive

    @staticmethod
    def process_archive(archive: MoomapArchive):
        for entry in archive.get_all():
            # Add original action info to assure same format for entries added from archive action
            info: dict[str, Any] = getattr(entry, "info", {})
            info["env_action"] = getattr(entry, "action", None)
            entry.info = info

    @staticmethod
    def parse_entry(entry: MoomapArchiveEntry, str_prefix: str = "") -> dict:
        x = entry.info["env_action"]
        y = entry.obs
        fit = entry._fitness
        bin = entry.bin_idx
        return {
            f"{str_prefix}x": x,
            f"{str_prefix}y": y,
            f"{str_prefix}fitness": fit,
            f"{str_prefix}bin": bin,
        }

    @staticmethod
    def run(config: MoomapXConfig, out_dir: str = "."):
        env_config = config.generate_env_config()
        config_dict = env_config.to_dict()  # Ensure the config is serializable
        env_list = []
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

            env_id = uuid.uuid4().hex
            env_info = {
                "id": env_id,
                "name": env.name,
                "func_id": env.func_id,
                "inst": env.inst,
                "ref_dirs": archive.binner.ref_dirs.tolist(),
                "ideal": archive.binner.ideal.tolist(),
                "nadir": archive.binner.nadir.tolist(),
            }
            env_list.append(env_info)
            res_list = []

            # Now measure the success of archive actions
            for _ in range(config.num_trials):
                action = (
                    waa_env.action_space.sample()
                )  # Sample from the archive action space
                res_dict = {f"a{i}_bin": v for i, v in enumerate(action)}
                res_dict["id"] = env_id
                # add prefix to all keys in archive_stats
                archive_stats = {
                    f"old_{k}": v for k, v in archive.get_archive_stats().items()
                }
                res_dict.update(archive_stats)
                try:
                    _obs, reward, _terminated, _truncated, info = waa_env.step(action)
                    for i, p in enumerate(info["parents"]):
                        res_dict.update(MoomapX.parse_entry(p, str_prefix=f"p{i}_"))
                    child = archive.last_entry
                    res_dict.update(MoomapX.parse_entry(child, str_prefix="c_"))
                    res_dict["reward"] = reward
                    res_dict["added"] = child.added
                    archive_stats = {
                        f"new_{k}": v for k, v in archive.get_archive_stats().items()
                    }
                    res_dict.update(archive_stats)
                except AssertionError as e:
                    assert str(e) == "Not enough parents found in the archive"
                # TODO: add HV as a metric to the archive stats
                res_list.append(res_dict)
            env.close()
            df = pd.DataFrame(res_list)
            df.to_csv(os.path.join(out_dir, f"results_{env_id}.csv"), index=False)
        # save experiment and envs info
        experiment_info = {
            "config": config_dict,
            "envs": env_list,
        }
        with open(os.path.join(out_dir, "experiment_info.json"), "w") as f:
            json.dump(
                experiment_info,
                f,
                indent=4,
                default=lambda x: x.tolist() if isinstance(x, np.ndarray) else x.item(),
            )
        return experiment_info


# TODO: Postprocessing that also determines distances
# TODO: add HV as metric to archive stats by implementing an unbounded HV archive + non_dominated
#
#
def parse_args():
    parser = argparse.ArgumentParser(description="Run MoomapX experiment")
    parser.add_argument(
        "--config_file", type=str, help="Path to the configuration file"
    )
    args = parser.parse_args()
    return args


def main(config_file: str, out_dir: str = "."):
    with open(config_file, "r") as f:
        run_config = json.load(f)
    config = CF.from_dict(run_config, context=globals())
    experiment_id = uuid.uuid4().hex
    # Make folder for experiment results for this experiment run
    out_dir = os.path.join(out_dir, f"x_{experiment_id}")
    os.makedirs(out_dir, exist_ok=False)
    MoomapX.run(config, out_dir=out_dir)
    return out_dir


if __name__ == "__main__":
    args = parse_args()
    main(args.config_file, out_dir=".")
