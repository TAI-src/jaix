import argparse
import copy
import json
import numbers
import os
import uuid
from enum import Enum

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
from cobi_config_generator import names as cobi_names


class Crossover(Enum):
    UNIFORM = "uniform"
    ONE_POINT = "one_point"
    TWO_POINT = "two_point"
    ARITHMETIC = "arithmetic"


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
        crossover: Crossover = Crossover.UNIFORM,
        num_parents: int = 2,
    ):
        super().__init__()
        self.num_independent_runs = num_independent_runs
        self.num_prefill_samples = num_prefill_samples
        self.num_offspring = num_offspring
        self.num_generations = num_generations
        self.seed = seed
        self.mode = mode
        self.num_parents = num_parents
        self.rng = np.random.default_rng(seed)
        self.independent_run_seeds = self.rng.integers(
            0, 2**32 - 1, size=num_independent_runs
        )
        self.crossover = crossover
        self.mo_archive_kwargs = (
            mo_archive_kwargs if mo_archive_kwargs is not None else {}
        )

    def update_defaults(self, problem: StaticProblem):
        self.mo_archive_config = NSGA3ExperimentConfig.create_mo_archive_config(
            problem, **self.mo_archive_kwargs
        )
        assert isinstance(self.mo_archive_config, MOArchiveConfig)
        assert isinstance(self.mo_archive_config.num_refpoints, int)
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
        # All non-constrained RE problems are included in the list. The constrained ones are excluded for now.
        problems = [REProblem(REProblemConfig(), i) for i in range(16)]
        return problems

    @staticmethod
    def cobi_problem_list():
        cobi_configs = get_cobi_configs()
        problems = [CobiProblem(config, inst=0) for config in cobi_configs]
        for i, problem in enumerate(problems):
            problem.name = cobi_names[i]

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
        assert (
            isinstance(config.num_prefill_samples, int)
            and config.num_prefill_samples > 0
        )
        for _ in range(config.num_prefill_samples):
            random_x = config.rng.uniform(problem.lower_bounds, problem.upper_bounds)
            _y_noise, y_raw = env.func(random_x)
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
            found_entry = archive.get(entry)
        assert found_entry is not None, "Entry not found in archive"
        assert isinstance(found_entry, MOEvalEntry)
        return {
            "x": found_entry.x,
            "y": found_entry.y,
            "rank": found_entry.rank,
            "sec_score": found_entry.secondary_score,
            "niche": found_entry.niche,
            "dist_to_ref": found_entry.dist_to_ref,
            "dist_to_ideal": found_entry.dist_to_ideal,
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
        assert len(parents) == config.num_parents
        offspring: np.ndarray | np.float32
        if config.crossover == Crossover.UNIFORM:
            if all(isinstance(p, numbers.Number) for p in parents):
                offspring = np.float32(np.mean(parents))
            elif all(isinstance(p, np.ndarray) for p in parents):
                offspring = np.asarray(np.mean(parents, axis=0), dtype=np.float32)
        else:
            raise NotImplementedError(
                f"Crossover method {config.crossover} not implemented"
            )
        return offspring

    @staticmethod
    def create_families(config: NSGA3ExperimentConfig, archive: MOArchive):
        families = []
        assert isinstance(config.num_offspring, int) and config.num_offspring > 0
        for _ in range(config.num_offspring):
            p_info = NSGA3Experiment.select_parents(config, archive)
            px_list = [p["x"] for p in p_info]

            offspring = NSGA3Experiment.create_offspring(px_list, config)
            fam_dict = {}
            for i, p in enumerate(p_info):
                fam_dict[f"parent_{i}"] = p
            fam_dict["offspring"] = offspring
            families.append(fam_dict)
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
            else:
                info_dict["added"] = False
            offspring_info.append(info_dict)

        return offspring_info, entries

    @staticmethod
    def run_single(
        o_config: NSGA3ExperimentConfig,
        out_dir: str = ".",
        problem_idx: list[int] | None = None,
    ):
        os.makedirs(out_dir, exist_ok=False)
        o_config.rng = np.random.default_rng(o_config.seed)
        files = []
        problem_list = o_config.generate_problem_list(o_config.mode)
        if problem_idx is not None:
            problem_list = [problem_list[i] for i in problem_idx]
        for problem in problem_list:
            config = copy.deepcopy(o_config)
            config.update_defaults(problem)
            results = []
            problem_dict = {"problem": str(problem), "seed": config.seed}
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
            df = pd.json_normalize(results, sep="_")
            file_name = f"{out_dir}/results_{problem!s}.csv"
            df.to_csv(file_name, index=False)
            files.append(file_name)
            # Also save the config to a json file
            config_dict = config.to_dict()
            file_name = f"{out_dir}/config_{problem!s}.json"
            with open(file_name, "w") as f:
                json.dump(
                    config_dict,
                    f,
                    indent=4,
                    default=lambda x: (
                        x.tolist() if isinstance(x, np.ndarray) else x.item()
                    ),
                )
            files.append(file_name)
        return files

    @staticmethod
    def run(
        config: NSGA3ExperimentConfig,
        out_dir: str = ".",
        problem_idx: list[int] | None = None,
    ):
        out_files = []
        exp_id = uuid.uuid4().hex
        for irun in config.independent_run_seeds:
            out_dir_run = f"{out_dir}/x_{exp_id}/r_{irun}"
            run_config = copy.deepcopy(config)
            run_config.seed = irun
            res = NSGA3Experiment.run_single(
                o_config=run_config, out_dir=out_dir_run, problem_idx=problem_idx
            )
            out_files.extend(res)
        config_dict = config.to_dict()
        config_dict["out_files"] = out_files
        config_dict["exp_id"] = exp_id
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


def parse_args():
    """
    num_independent_runs: int,
    num_generations: int,
    num_prefill_samples: int = -1,
    num_offspring: int = -1,
    mode: str = "",
    seed: int | None = None,
    mo_archive_kwargs: dict | None = None,
    crossover: Crossover = Crossover.UNIFORM,
    num_parents: int = 2,
    """
    parser = argparse.ArgumentParser(description="Run NSGA3 experiment")
    parser.add_argument(
        "--num_independent_runs", type=int, help="Number of independent runs"
    )
    parser.add_argument("--num_generations", type=int, help="Number of generations")
    parser.add_argument(
        "--num_prefill_samples",
        type=int,
        default=-1,
        help="Number of prefill samples",
    )
    parser.add_argument(
        "--num_offspring", type=int, default=-1, help="Number of offspring"
    )
    parser.add_argument("--mode", type=str, default="", help="Mode: cobi or reproblem")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument(
        "--mo_archive_kwargs",
        type=json.loads,
        default={},
        help="MO archive kwargs as JSON string",
    )
    parser.add_argument(
        "--crossover",
        type=str,
        default=Crossover.UNIFORM.value,
        choices=[c.value for c in Crossover],
        help="Crossover method",
    )
    parser.add_argument(
        "--out_dir", type=str, default=".", help="Output directory for results"
    )
    parser.add_argument("--num_parents", type=int, default=2, help="Number of parents")
    args = parser.parse_args()
    return args


def main(args):
    config = NSGA3ExperimentConfig(
        num_independent_runs=args.num_independent_runs,
        num_generations=args.num_generations,
        num_prefill_samples=args.num_prefill_samples,
        num_offspring=args.num_offspring,
        mode=args.mode,
        seed=args.seed,
        mo_archive_kwargs=args.mo_archive_kwargs,
        crossover=Crossover(args.crossover),
        num_parents=args.num_parents,
    )
    NSGA3Experiment.run(config, out_dir=args.out_dir)


if __name__ == "__main__":
    args = parse_args()
    main(args)
