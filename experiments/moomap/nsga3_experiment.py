import copy
import json
import numbers
import os
import uuid

import numpy as np
import pandas as pd
from jaix.env.singular.ec_env import ECEnvironment, ECEnvironmentConfig
from jaix.env.utils.archive.entry_scorer import (
    EntryScorer,
    ReferenceVectorDistanceScorer,
)
from jaix.env.utils.archive.mo_archive import (
    KeepDominated,
    MOArchive,
    MOArchiveConfig,
    MOArchiveEntry,
)
from jaix.env.utils.mo_sizing import get_num_refpoints
from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from jaix.env.utils.problem.static_problem import StaticProblem
from ttex.config import Config

from cobi_config_generator import get_configs as get_cobi_configs


class NSGA3ExperimentConfig(Config):
    def __init__(
        self,
        num_independent_runs: int,
        num_generations: int,
        num_prefill_samples: int = -1,
        num_offspring: int = -1,
        mode: str = "",
        seed: int | None = None,
        mo_archive_kwargs: dict | None = None,
    ):
        super().__init__()
        self.num_independent_runs = num_independent_runs
        self.num_prefill_samples = num_prefill_samples
        self.num_offspring = num_offspring
        self.num_generations = num_generations
        self.seed = seed
        self.mode = mode
        self.num_parents = 2  # Number of parents for the uniform crossover action space
        self.rng = np.random.default_rng(seed)
        self.independent_run_seeds = self.rng.integers(
            0, 2**32 - 1, size=num_independent_runs
        )
        self.mo_archive_kwargs = (
            mo_archive_kwargs if mo_archive_kwargs is not None else {}
        )

    def update_defaults(self, problem: StaticProblem):
        self.mo_archive_config = NSGA3ExperimentConfig.create_mo_archive_config(
            problem, **self.mo_archive_kwargs
        )
        self.num_prefill_samples = (
            self.mo_archive_config.num_refpoints
            if self.num_prefill_samples < 0
            else self.num_prefill_samples
        )
        self.num_offspring = (
            self.mo_archive_config.num_refpoints
            if self.num_offspring < 0
            else self.num_offspring
        )

    @staticmethod
    def create_mo_archive_config(
        problem: StaticProblem,
        secondary_criterion_class: type[EntryScorer] = ReferenceVectorDistanceScorer,
        max_size: int | None = None,
        keep_dominated: KeepDominated = KeepDominated.ALL,
        only_new_entries: bool = False,
        hv_approx_samples: int | None = 262_144,
        num_refpoints: int | str = "original",
    ) -> MOArchiveConfig:
        if isinstance(num_refpoints, str):
            n_refpoints: int = get_num_refpoints(problem.num_objectives, num_refpoints)
        else:
            n_refpoints = num_refpoints
        max_size = n_refpoints if max_size is None else max_size

        config = MOArchiveConfig(
            archive_entry_class=MOEvalEntry,
            secondary_criterion_class=secondary_criterion_class,
            max_size=max_size,
            keep_dominated=keep_dominated,
            only_new_entries=only_new_entries,
            hv_approx_samples=hv_approx_samples,
            num_refpoints=n_refpoints,
        )
        return config

    @staticmethod
    def re_problem_list():
        problems = [REProblem(REProblemConfig(), i) for i in range(24)]
        return problems

    @staticmethod
    def cobi_problem_list():
        cobi_configs = get_cobi_configs()
        problems = [CobiProblem(config, inst=0) for config in cobi_configs]
        return problems

    @staticmethod
    def generate_problem_list(mode: str = ""):
        if mode == "cobi":
            return NSGA3ExperimentConfig.cobi_problem_list()
        elif mode == "reproblem":
            return NSGA3ExperimentConfig.re_problem_list()
        else:
            return (
                NSGA3ExperimentConfig.cobi_problem_list()
                + NSGA3ExperimentConfig.re_problem_list()
            )


class MOEvalEntry(MOArchiveEntry):
    def __init__(self, x: np.ndarray, y: np.ndarray):
        self.x = x
        self.y = y

    def parse(self) -> np.ndarray:
        return self.y


class NSGA3Experiment:

    @staticmethod
    def prefill_archive(problem: StaticProblem, config: NSGA3ExperimentConfig):
        env = ECEnvironment(ECEnvironmentConfig(budget_multiplier=1), func=problem)
        archive = MOArchive(
            config.mo_archive_config, env=env
        )  # env is not used in prefill
        entries = []
        for _ in range(config.num_prefill_samples):
            random_x = config.rng.uniform(problem.lower_bounds, problem.upper_bounds)
            _y_noise, y_raw = env(random_x)
            entry = MOEvalEntry(random_x, np.asarray(y_raw))
            entries.append(entry)
        archive.add(entries)
        return archive, entries

    @staticmethod
    def get_unbounded_archive(problem: StaticProblem, config: NSGA3ExperimentConfig):
        env = ECEnvironment(ECEnvironmentConfig(budget_multiplier=1), func=problem)
        archive_config = copy.deepcopy(config.mo_archive_config)
        archive_config.keep_dominated = (
            KeepDominated.NONE
        )  # only keep non-dominated entries for unbounded archive
        archive_config.max_size = None  # no size limit for unbounded archive
        archive_config.only_new_entries = False  # keep elites
        archive = MOArchive(archive_config, env=env)
        return archive

    @staticmethod
    def entry_info(archive: MOArchive, entry: MOEvalEntry | int):
        if isinstance(entry, MOEvalEntry):
            found_entry = entry
        else:
            found_entry = archive.get_entry(entry)
        return {
            "x": found_entry.x,
            "y": found_entry.y,
            "rank": found_entry.rank,
            "sec_score": found_entry.secondary_score,
            "niche": found_entry.niche,
            "dist_to_ref": found_entry.dist_to_ref,
        }

    @staticmethod
    def select_parents(config: NSGA3ExperimentConfig, archive: MOArchive):
        p_info = []
        # identify parents from the archive
        for _ in range(config.num_parents):
            parent_idx = int(config.rng.integers(0, archive.size))
            p_dict = NSGA3Experiment.entry_info(archive, parent_idx)
            p_info.append(p_dict)
        return p_info

    @staticmethod
    def create_offspring(
        parents: list[np.ndarray], config: NSGA3ExperimentConfig
    ) -> np.ndarray | np.float32:
        # TODO: Different crossover methods can be implemented here. For now, we will use uniform crossover.
        if all(isinstance(p, numbers.Number) for p in parents):
            offspring = np.float32(np.mean(parents))
        elif all(isinstance(p, np.ndarray) for p in parents):
            offspring = np.asarray(np.mean(parents, axis=0), dtype=np.float32)
        return offspring

    @staticmethod
    def create_families(config: NSGA3ExperimentConfig, archive: MOArchive):
        families = []
        for _ in range(config.num_offspring):
            p_info = NSGA3Experiment.select_parents(config, archive)
            px_list = [p["x"] for p in p_info]

            offspring = NSGA3Experiment.create_offspring(px_list, config)
            families.append({"parents": p_info, "offspring": offspring})
        return families

    @staticmethod
    def try_add_offspring(
        problem: StaticProblem, archive: MOArchive, offspring: list[np.ndarray]
    ):
        offspring_info = []
        entries = []
        for off in offspring:
            _y_noise, y_raw = problem(off)
            entry = MOEvalEntry(off, np.asarray(y_raw))
            entries.append(entry)
        archive.add(entries)
        for entry in entries:
            info_dict = NSGA3Experiment.entry_info(archive, entry)
            if entry in archive.archived_entries:
                info_dict["added"] = True
            offspring_info.append(info_dict)

        return offspring_info, entries

    @staticmethod
    def run_single(config: NSGA3ExperimentConfig, out_dir: str = "."):
        for problem in config.generate_problem_list():
            config.update_defaults(problem)
            results = []
            problem_dict = {"problem": str(problem)}
            archive, all_evaluated = NSGA3Experiment.prefill_archive(problem, config)
            unbounded_archive = NSGA3Experiment.get_unbounded_archive(problem, config)
            unbounded_archive.add(all_evaluated)

            for gen in range(config.num_generations):
                prev_archive_stats = archive.get_archive_stats()
                prev_archive_stats["unbounded_hv"] = unbounded_archive.score
                families = NSGA3Experiment.create_families(config, archive)
                offspring = [fam["offspring"] for fam in families]
                offspring_info, all_offspring = NSGA3Experiment.try_add_offspring(
                    problem, archive, offspring
                )
                unbounded_archive.add(all_offspring)
                new_archive_stats = archive.get_archive_stats()
                new_archive_stats["unbounded_hv"] = unbounded_archive.score
                gen_dict = {
                    "archive_stats_before": prev_archive_stats,
                    "archive_stats_after": new_archive_stats,
                    "generation": gen,
                }
                gen_dict.update(problem_dict)
                for fam_info, off_info in zip(families, offspring_info):
                    fam_info["offspring"] = off_info
                    fam_info.update(gen_dict)
                    results.append(fam_info)
            # create data frame and save to csv
            df = pd.DataFrame(results)
            file_name = f"{out_dir}/results_{problem!s}.csv"
            df.to_csv(file_name, index=False)
            # Also save the config to a json file
            config_dict = config.to_dict()
            with open(f"{out_dir}/config_{problem!s}.json", "w") as f:
                json.dump(
                    config_dict,
                    f,
                    indent=4,
                    default=lambda x: (
                        x.tolist() if isinstance(x, np.ndarray) else x.item()
                    ),
                )
        return out_dir

    @staticmethod
    def run(config: NSGA3ExperimentConfig, out_dir: str = "."):
        out_dirs = []
        exp_id = uuid.uuid4().hex
        for irun in config.independent_run_seeds:
            out_dir_run = f"{out_dir}/x_{exp_id}/r_{irun}"
            os.makedirs(out_dir_run, exist_ok=True)
            run_config = copy.deepcopy(config)
            run_config.rng = np.random.default_rng(irun)
            res = NSGA3Experiment.run_single(config=config, out_dir=out_dir_run)
            out_dirs.append(res)
        config_dict = config.to_dict()
        config_dict["out_dirs"] = out_dirs
        # Print config to a file in the out_dir
        file_name = f"{out_dir}/x_{exp_id}/config.json"
        with open(file_name, "w") as f:
            json.dump(
                config_dict,
                f,
                indent=4,
                default=lambda x: x.tolist() if isinstance(x, np.ndarray) else x.item(),
            )
        return config_dict
