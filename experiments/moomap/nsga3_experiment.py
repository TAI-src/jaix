import copy
import json
import os
import uuid

import numpy as np
import pandas as pd
from jaix.env.singular.ec_env import ECEnvironment, ECEnvironmentConfig
from jaix.env.utils.archive.mo_archive import (
    KeepDominated,
    MOArchive,
    MOArchiveEntry,
)
from jaix.env.utils.problem.static_problem import StaticProblem

from config_nsga3x import (
    Crossover,
    MOEvalEntry,
    NSGA3ExperimentConfig,
    generate_config,
    parse_args,
)
from utils_problems import generate_problem_list


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
        assert isinstance(found_entry, MOEvalEntry) and isinstance(
            found_entry, MOArchiveEntry
        )
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
        parents: list[np.ndarray],
        config: NSGA3ExperimentConfig,
    ) -> np.ndarray | np.float32:
        assert len(parents) == config.num_parents
        if config.crossover == Crossover.UNIFORM:
            weight_vec = config.rng.uniform(0, 1, size=len(parents))
            normalized_weights = weight_vec / np.sum(weight_vec)
            offspring = np.sum(
                [w * p for w, p in zip(normalized_weights, parents)], axis=0
            )
            return offspring
        raise NotImplementedError(
            f"Crossover method {config.crossover} not implemented"
        )

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
        problem_list = generate_problem_list(problem_idx)
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
                if gen % 100 == 0 or gen == config.num_generations - 1:
                    print(
                        f"Problem: {problem}, Generation: {gen}, Archive size: {archive.size}, results collected: {len(results)}"
                    )
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


def main(args):
    config = generate_config(args)
    NSGA3Experiment.run(config, out_dir=args.out_dir, problem_idx=args.problem_idx)


if __name__ == "__main__":
    args = parse_args()
    main(args)
