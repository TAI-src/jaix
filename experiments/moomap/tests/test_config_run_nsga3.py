import pytest

from config_run_nsga3 import get_batches, parse_args


def test_get_batches_full():
    settings, num_batches = get_batches(
        batch_id=None, n_runs=2, problem_ids=None, seed=42
    )
    assert len(settings) == num_batches
    assert num_batches == 2 * 23 * 2
    unique_seeds = {setting["seed"] for setting in settings}
    assert len(unique_seeds) == 2

    settings2, nb2 = get_batches(batch_id=[3], n_runs=2, problem_ids=None, seed=42)
    assert len(settings2) == 1
    assert settings2[0]["id"] == 3
    assert str(settings2[0]["problem"]) == str(settings[3]["problem"])
    assert settings2[0]["seed"] == settings[3]["seed"]
    assert settings2[0]["static_ref"] == settings[3]["static_ref"]
    assert nb2 == num_batches


def test_with_static_ref_only():
    settings, num_batches = get_batches(
        batch_id=None, n_runs=2, problem_ids=None, static_ref=True, seed=42
    )
    assert all(setting["static_ref"] is True for setting in settings)
    assert len(settings) == 2 * 23 * 1
    assert num_batches == 2 * 23 * 1


test_scenarios = [
    (["--static_ref"], True),
    ([], None),
    (["--not_static_ref"], False),
    (["--static_ref", "--not_static_ref"], None),
]


@pytest.mark.parametrize("args, expected_static_ref", test_scenarios)
def test_parse_args(monkeypatch, args, expected_static_ref):
    test_args = [
        "script_name",
        "--n_gen",
        "10",
        "--seed",
        "42",
        "--problem_ids",
        "0",
        "1",
    ]
    test_args.extend(args)
    monkeypatch.setattr("sys.argv", test_args)
    args = parse_args()
    assert args.n_gen == 10
    assert args.seed == 42
    assert args.static_ref is expected_static_ref
    assert args.problem_ids == [0, 1]
