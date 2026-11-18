from mating_kernel.exp.utils import batch
import pytest
import argparse
from mating_kernel.exp.utils.parse_args import (
    batched_experiment_parser,
    experiment_parser,
    exp_parser,
    results_parser,
)


def test_experiment_parser_defaults(tmp_path):
    parser = batched_experiment_parser()

    args = parser.parse_args(["--out_dir", str(tmp_path / "results")])

    assert args.cobi is False
    assert args.re is False
    assert args.num_objectives is None
    assert args.constrained is False

    assert args.seed is None
    assert args.reps == 30
    assert args.num_batches == 1
    assert args.out_dir == str(tmp_path / "results")
    assert args.nth is None


def test_experiment_parser_parses_all_arguments(tmp_path):
    parser = batched_experiment_parser()

    args = parser.parse_args(
        [
            "--cobi",
            "--re",
            "--num_objectives",
            "2",
            "3",
            "--constrained",
            "--seed",
            "123",
            "--reps",
            "10",
            "--num_batches",
            "4",
            "--out_dir",
            str(tmp_path / "results"),
            "--nth",
            "0",
            "3",
            "--group_by",
            "pid",
            "sid",
        ]
    )

    assert args.cobi is True
    assert args.re is True
    assert args.num_objectives == [2, 3]
    assert args.constrained is True

    assert args.seed == 123
    assert args.reps == 10
    assert args.num_batches == 4
    assert args.out_dir == str(tmp_path / "results")
    assert args.nth == [0, 3]
    assert args.group_by == ["pid", "sid"]
    assert args.mode == "run"
    assert args.skip_existing is False
    assert args.force_recompute is False


def test_exp_parser_with_options(tmp_path):
    parser = batched_experiment_parser()

    args = parser.parse_args(
        [
            "--out_dir",
            str(tmp_path / "results"),
            "--mode",
            "check",
            "--skip_existing",
            "--force_recompute",
        ]
    )

    assert args.out_dir == str(tmp_path / "results")
    assert args.mode == "check"
    assert args.skip_existing is True
    assert args.force_recompute is True


def test_experiment_parser_requires_out_dir():
    parser = batched_experiment_parser()

    with pytest.raises(SystemExit):
        parser.parse_args([])


def test_experiment_parser_rejects_unknown_arguments(tmp_path):
    parser = batched_experiment_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "--out_dir",
                str(tmp_path / "results"),
                "results",
                "--does-not-exist",
            ]
        )


def test_results_parser(tmp_path):
    parser = results_parser()

    args = parser.parse_args(
        [
            "--input_dirs",
            str(tmp_path / "input1"),
            "--out_dir",
            str(tmp_path / "results"),
            "--mode",
            "results",
        ]
    )
    assert args.input_dirs == [str(tmp_path / "input1")]
    assert args.out_dir == str(tmp_path / "results")
    assert args.mode == "results"


@pytest.mark.parametrize("mode", [None, "check", "pp", "run", "results"])
def test_experiments_parser(mode):
    parser = experiment_parser(mode=mode)
    assert isinstance(parser, argparse.ArgumentParser)
    args, unknown = parser.parse_known_args(
        ["--mode", mode if mode is not None else "run", "--out_dir", "results"]
    )
    assert args.mode == (mode if mode is not None else "run")
