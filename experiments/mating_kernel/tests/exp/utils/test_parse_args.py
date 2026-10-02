import pytest

from mating_kernel.exp.utils.parse_args import experiment_parser


def test_experiment_parser_defaults():
    parser = experiment_parser()

    args = parser.parse_args(
        [
            "--out_dir",
            "results",
        ]
    )

    assert args.cobi is False
    assert args.re is False
    assert args.num_objectives is None
    assert args.constrained is False

    assert args.seed is None
    assert args.reps == 30
    assert args.num_batches == 1
    assert args.out_dir == "results"
    assert args.nth is None


def test_experiment_parser_parses_all_arguments():
    parser = experiment_parser()

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
            "results",
            "--nth",
            "0",
            "3",
        ]
    )

    assert args.cobi is True
    assert args.re is True
    assert args.num_objectives == [2, 3]
    assert args.constrained is True

    assert args.seed == 123
    assert args.reps == 10
    assert args.num_batches == 4
    assert args.out_dir == "results"
    assert args.nth == [0, 3]


def test_experiment_parser_requires_out_dir():
    parser = experiment_parser()

    with pytest.raises(SystemExit):
        parser.parse_args([])


def test_experiment_parser_rejects_unknown_arguments():
    parser = experiment_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "--out_dir",
                "results",
                "--does-not-exist",
            ]
        )
