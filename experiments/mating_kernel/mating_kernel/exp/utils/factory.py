import argparse
from enum import Enum
from pathlib import Path
from typing import Any

from ttex.config import Config

from mating_kernel.problems.mk_suite import MKSuiteConfig


class ExperimentMode(Enum):
    RUN = "run"
    CHECK = "check"
    PP = "pp"
    RESULTS = "results"


class ExperimentConfig(Config):
    def __init__(
        self,
        mode: ExperimentMode,
        out_dir: str | Path = Path("."),
        settings: dict[str, Any] | None = None,
    ):
        super().__init__()
        self.mode = mode
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.settings = settings if settings is not None else {}


class ResultsConfig(ExperimentConfig):
    def __init__(
        self,
        input_dirs: list[str] | list[Path],
        out_dir: str | Path,
        settings: dict[str, Any] | None = None,
    ):
        super().__init__(ExperimentMode.RESULTS, out_dir=out_dir, settings=settings)
        self.input_dirs = [Path(d) for d in input_dirs]


class BatchedExperimentConfig(ExperimentConfig):
    def __init__(
        self,
        suite_config: MKSuiteConfig,
        settings: dict[str, list],
        reps: int = 1,
        num_batches: int | None = None,
        seed: int | None = None,
        out_dir: str | Path = Path("."),
        mode: ExperimentMode = ExperimentMode.RUN,
        group_by: list[str] | None = None,
        skip_existing: bool = False,
        force_recompute: bool = False,
    ):
        super().__init__(mode=mode, out_dir=out_dir, settings=settings)
        self.suite_config = suite_config
        self.reps = reps
        self.num_batches = num_batches
        self.seed = seed
        self.group_by = group_by
        self.skip_existing = skip_existing
        self.force_recompute = force_recompute


def parse_mk_suite_args(args: argparse.Namespace) -> MKSuiteConfig:
    return MKSuiteConfig(
        cobi=args.cobi,
        re=args.re,
        num_objectives=args.num_objectives,
        constrained=args.constrained,
    )


def parse_config(
    mode: ExperimentMode, args: argparse.Namespace, settings: dict[str, Any]
) -> ResultsConfig | BatchedExperimentConfig:
    if mode == ExperimentMode.RESULTS:
        return parse_results_config(args, settings)
    else:
        return parse_experiment_config(args, settings)


def parse_results_config(
    args: argparse.Namespace, settings: dict[str, Any]
) -> ResultsConfig:
    return ResultsConfig(
        input_dirs=args.input_dirs,
        out_dir=args.out_dir,
        settings=settings,
    )


def parse_experiment_config(
    args: argparse.Namespace, settings: dict[str, Any]
) -> BatchedExperimentConfig:
    suite_config = parse_mk_suite_args(args)
    return BatchedExperimentConfig(
        suite_config=suite_config,
        settings=settings,
        reps=args.reps,
        num_batches=args.num_batches,
        seed=args.seed,
        out_dir=args.out_dir,
        mode=ExperimentMode(args.mode),
        group_by=args.group_by,
        skip_existing=args.skip_existing,
        force_recompute=args.force_recompute,
    )
