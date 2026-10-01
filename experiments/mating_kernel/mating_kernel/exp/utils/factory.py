from typing import Any
from dataclasses import dataclass
import argparse
from mating_kernel.problems.mk_suite import MKSuiteConfig
from mating_kernel.exp.utils.parse_args import mk_suite_parser, experiment_parser
from mating_kernel.experiments.experiment import ExperimentConfig


@dataclass
class SettingArg:
    args: tuple[str, ...]
    kwargs: dict[str, Any]


def convert_settings(
    settings: list[SettingArg],
) -> tuple[argparse.ArgumentParser, list[str]]:
    parser = argparse.ArgumentParser(
        add_help=False, parents=[mk_suite_parser(), experiment_parser()]
    )
    setting_names = []
    for setting in settings:
        action = parser.add_argument(*setting.args, **setting.kwargs)
        setting_names.append(action.dest)
    return parser, setting_names


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
