import argparse
import json
from enum import Enum

import numpy as np
from jaix.env.utils.archive.entry_scorer import (
    EntryScorer,
    ReferenceVectorDistanceScorer,
)
from jaix.env.utils.archive.mo_archive import (
    KeepDominated,
    MOArchiveConfig,
    MOArchiveEntry,
)
from jaix.env.utils.mo_sizing import get_num_refpoints
from jaix.env.utils.problem.static_problem import StaticProblem
from ttex.config import Config


class MOEvalEntry(MOArchiveEntry):
    def __init__(self, x: np.ndarray, y: np.ndarray):
        self.x = x
        self.y = y

    def parse(self) -> np.ndarray:
        return self.y


class Crossover(Enum):
    UNIFORM = "uniform"
    ONE_POINT = "one_point"
    TWO_POINT = "two_point"
    ARITHMETIC = "arithmetic"


def parse_args():

    parser = argparse.ArgumentParser(description="Run NSGA3 experiment")
    parser.add_argument(
        "--num_independent_runs",
        type=int,
        help="Number of independent runs",
        required=True,
    )
    parser.add_argument(
        "--num_generations", type=int, help="Number of generations", required=True
    )
    parser.add_argument(
        "--num_prefill_samples",
        type=int,
        default=-1,
        help="Number of prefill samples",
    )
    parser.add_argument(
        "--num_offspring", type=int, default=-1, help="Number of offspring"
    )
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
    parser.add_argument(
        "--problem_idx", type=int, nargs="*", default=None, help="Problem indices"
    )
    parser.add_argument("--num_parents", type=int, default=2, help="Number of parents")
    args = parser.parse_args()
    return args


class NSGA3ExperimentConfig(Config):
    def __init__(
        self,
        num_independent_runs: int,
        num_generations: int,
        num_prefill_samples: int = -1,
        num_offspring: int = -1,
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
        self.num_parents = num_parents
        assert num_parents > 0, "Number of parents must be greater than 0"
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


def generate_config(args):
    config = NSGA3ExperimentConfig(
        num_independent_runs=args.num_independent_runs,
        num_generations=args.num_generations,
        num_prefill_samples=args.num_prefill_samples,
        num_offspring=args.num_offspring,
        seed=args.seed,
        mo_archive_kwargs=args.mo_archive_kwargs,
        crossover=Crossover(args.crossover),
        num_parents=args.num_parents,
    )
    return config
