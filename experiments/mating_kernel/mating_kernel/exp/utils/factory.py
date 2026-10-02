import argparse
from pathlib import Path
from typing import Any

from ttex.config import Config

from mating_kernel.problems.mk_suite import MKSuiteConfig


class ExperimentConfig(Config):
    def __init__(
        self,
        suite_config: MKSuiteConfig,
        settings: dict[str, list],
        reps: int = 1,
        num_batches: int | None = None,
        seed: int | None = None,
        out_dir: str | Path = Path("."),
    ):
        super().__init__()
        self.suite_config = suite_config
        self.reps = reps
        self.settings = settings
        self.num_batches = num_batches
        self.seed = seed
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)


def parse_mk_suite_args(args: argparse.Namespace) -> MKSuiteConfig:
    return MKSuiteConfig(
        cobi=args.cobi,
        re=args.re,
        num_objectives=args.num_objectives,
        constrained=args.constrained,
    )


def parse_experiment_config(
    args: argparse.Namespace, settings: dict[str, Any]
) -> ExperimentConfig:
    suite_config = parse_mk_suite_args(args)
    return ExperimentConfig(
        suite_config=suite_config,
        settings=settings,
        reps=args.reps,
        num_batches=args.num_batches,
        seed=args.seed,
        out_dir=args.out_dir,
    )
