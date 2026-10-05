import argparse
from pathlib import Path

from mating_kernel.exp.utils.factory import (
    ExperimentConfig,
    parse_experiment_config,
    parse_mk_suite_args,
    ExperimentMode,
)
from mating_kernel.problems.mk_suite import MKSuiteConfig


def test_mk_suite_parser():
    args = argparse.Namespace(
        cobi=True,
        re=False,
        num_objectives=[2, 3],
        constrained=True,
    )
    mk_suite_config = parse_mk_suite_args(args)
    assert mk_suite_config.cobi is True
    assert mk_suite_config.re is False
    assert mk_suite_config.num_objectives == [2, 3]
    assert mk_suite_config.constrained is True
    assert isinstance(mk_suite_config, MKSuiteConfig)


def test_experiment_config_parser():
    args = argparse.Namespace(
        cobi=True,
        re=True,
        num_objectives=[2],
        constrained=False,
        reps=10,
        num_batches=5,
        seed=42,
        out_dir="results",
        mode="run",
        group_by=["pid", "sid"],
        skip_existing=True,
        force_recompute=True,
    )
    settings = {"some_setting": "value"}
    experiment_config = parse_experiment_config(args, settings)
    assert experiment_config.suite_config.cobi is True
    assert experiment_config.suite_config.re is True
    assert experiment_config.suite_config.num_objectives == [2]
    assert experiment_config.suite_config.constrained is False
    assert experiment_config.reps == 10
    assert experiment_config.num_batches == 5
    assert experiment_config.seed == 42
    assert experiment_config.out_dir == Path("results")
    assert experiment_config.settings == settings
    assert experiment_config.mode == ExperimentMode.RUN
    assert experiment_config.group_by == ["pid", "sid"]
    assert experiment_config.skip_existing is True
    assert experiment_config.force_recompute is True
    assert isinstance(experiment_config, ExperimentConfig)
