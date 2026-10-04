import shutil

import pytest
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)

from mating_kernel.exp.utils.batch import Batch
from mating_kernel.exp.utils.batcher import Batcher
from mating_kernel.problems.problem_info import ProblemInfo


def make_batch(tmp_path, **kwargs) -> Batch:
    problem = REProblem(REProblemConfig(), 0)
    pinfo = ProblemInfo(problem)
    defaults = {
        "pid": "p1",
        "sid": 0,
        "rep": 0,
        "parent_dir": tmp_path,
        "problem": problem,
        "pinfo": pinfo,
    }
    defaults.update(kwargs)
    return Batch(**defaults)


def sim_run(batch: Batch):
    # Simulate running the batch by creating a result file in the output directory
    batch.out_dir.mkdir(parents=True, exist_ok=True)
    with open(batch.out_dir / "result.txt", "w") as f:
        f.write("done")


@pytest.mark.parametrize("skip_existing", [True, False])
def test_seed_batches(skip_existing, tmp_path):
    batch = make_batch(tmp_path)

    result = Batcher.seed_batches([batch], reps=3, seed=123, skip_existing=False)

    assert len(result) == 3
    assert [b.rep for b in result] == [0, 1, 2]
    assert all(b.seed is not None for b in result)
    # simulated run experiments by creating files in outdir
    for b in result:
        sim_run(b)
    # remove out_dir for batch2 and see that it is picked up again
    seeds = [b.seed for b in result]
    shutil.rmtree(result[1].out_dir)
    result1 = Batcher.seed_batches(
        [batch], reps=3, seed=123, skip_existing=skip_existing
    )
    if skip_existing:
        assert len(result1) == 1
        assert [b.rep for b in result1] == [1]
        assert result1[0].seed == seeds[1]
    else:
        assert len(result1) == 3
        assert [b.rep for b in result1] == [0, 1, 2]
        assert [b.seed for b in result1] == seeds

    result2 = Batcher.seed_batches(
        [batch],
        reps=3,
        seed=124,
    )
    assert len(result2) == 3
    assert [b.rep for b in result2] == [0, 1, 2]
    assert all(b.seed is not None for b in result2)
    assert [b.seed for b in result2] != [b.seed for b in result]


def test_multiple_batches(tmp_path):
    batch1 = make_batch(tmp_path, sid=0, test="a")
    batch2 = make_batch(tmp_path, sid=1, test="b")
    batches = [batch1, batch2]

    result = Batcher.seed_batches(
        batches,
        reps=2,
        seed=123,
    )

    assert len(result) == 4
    assert [b.sid for b in result] == [0, 1, 0, 1]
    assert [b.rep for b in result] == [0, 0, 1, 1]

    # simulate runnign first batch
    sim_run(result[0])
    result2 = Batcher.seed_batches(
        batches,
        reps=2,
        seed=123,
    )
    assert len(result2) == 3
    assert [b.sid for b in result2] == [0, 1, 1]
    assert [b.rep for b in result2] == [1, 0, 1]


def test_split_batches_none():
    batches = [1, 2, 3]

    result = Batcher.split_batches(batches)

    assert result == [[1], [2], [3]]


def test_split_batches_one():
    batches = [1, 2, 3, 4]

    result = Batcher.split_batches(batches, num_batches=1)

    assert result == [[1, 2, 3, 4]]


def test_split_batches_even():
    batches = list(range(6))

    result = Batcher.split_batches(batches, num_batches=3)

    assert result == [
        [0, 1],
        [2, 3],
        [4, 5],
    ]


def test_split_batches_uneven():
    batches = list(range(10))

    result = Batcher.split_batches(batches, num_batches=3)

    assert result == [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [8, 9],
    ]


def test_split_batches_more_batches_than_items():
    batches = [1, 2, 3]

    result = Batcher.split_batches(batches, num_batches=10)

    assert result == [
        [1],
        [2],
        [3],
    ]


def test_split_batches_rejects_zero():
    with pytest.raises(AssertionError):
        Batcher.split_batches([1, 2, 3], num_batches=0)


def test_create_combinations(tmp_path):
    from mating_kernel.problems.mk_suite import MKSuite, MKSuiteConfig

    mk_config = MKSuiteConfig(cobi=True, re=True, constrained=False, num_objectives=[4])
    mk_suite = MKSuite(mk_config)
    assert len(mk_suite.problems) == 2  # There should be 2 problems with 4 objectives
    settings = {"setting1": [1, 2], "setting2": ["a", "b", "c"]}
    combinations = Batcher.create_combinations(mk_config, settings, tmp_path)
    assert len(combinations) == len(mk_suite.problems) * len(
        settings["setting1"]
    ) * len(settings["setting2"])
    fbatch = combinations[0]
    assert fbatch.pid == "RE41"
    assert fbatch.sid == 0
    assert fbatch.settings["setting1"] == 1
    assert fbatch.settings["setting2"] == "a"
    assert fbatch.problem.num_objectives == 4
    assert fbatch.pinfo.num_objectives == 4
    assert fbatch.rep is None
    assert fbatch.seed is None


def test_group_batches(tmp_path):
    batch1 = make_batch(tmp_path, sid=0, test="a")
    batch2 = make_batch(tmp_path, sid=1, test="b")
    batch3 = make_batch(tmp_path, sid=0, test="c")
    batches = [batch1, batch2, batch3]

    grouped = Batcher.group_batches(batches, group_by=["sid"])

    assert len(grouped) == 2
    assert len(grouped["0"]) == 2
    assert len(grouped["1"]) == 1
    assert grouped["0"][0].test == "a"
    assert grouped["0"][1].test == "c"
    assert grouped["1"][0].test == "b"


def test_group_batches_multiple_keys(tmp_path):
    batch1 = make_batch(tmp_path, sid=0, test="a")
    batch2 = make_batch(tmp_path, sid=1, test="b")
    batch3 = make_batch(tmp_path, sid=0, test="c")
    batch4 = make_batch(tmp_path, sid=1, test="d")
    batches = [batch1, batch2, batch3, batch4]

    grouped = Batcher.group_batches(batches, group_by=["sid", "test"])

    assert len(grouped) == 4
    assert grouped["0_a"][0].test == "a"
    assert grouped["0_c"][0].test == "c"
    assert grouped["1_b"][0].test == "b"
    assert grouped["1_d"][0].test == "d"


def test_create_batches(tmp_path):
    from mating_kernel.problems.mk_suite import MKSuiteConfig

    mk_config = MKSuiteConfig(cobi=True, re=True, constrained=False, num_objectives=[4])
    settings = {"setting1": [1, 2], "setting2": ["a", "b", "c"]}
    batches = Batcher.create_batches(
        suite_config=mk_config,
        settings=settings,
        reps=2,
        num_batches=3,
        seed=123,
        exp_dir=tmp_path,
        filter_bids=None,
        skip_existing=False,
    )
    # There should be 2 problems * 2 * 3 = 12 combinations, each with 2 reps = 24 batches
    assert len(batches) == 3  # Split into 3 batches
    total_batches = sum(len(batch_group) for batch_group in batches)
    assert total_batches == 24


def test_create_batches_with_filter(tmp_path):
    from mating_kernel.problems.mk_suite import MKSuiteConfig

    mk_config = MKSuiteConfig(cobi=True, re=True, constrained=False, num_objectives=[4])
    settings = {"setting1": [1, 2], "setting2": ["a", "b", "c"]}
    batches = Batcher.create_batches(
        suite_config=mk_config,
        settings=settings,
        reps=2,
        num_batches=3,
        seed=123,
        exp_dir=tmp_path,
        filter_bids=[0, 2],  # Only take the first and third batch groups
        skip_existing=False,
    )
    assert len(batches) == 2  # Only two batch groups should be returned
    total_batches = sum(len(batch_group) for batch_group in batches)
    assert (
        total_batches == 16
    )  # Each group should have 8 batches (24 total / 3 groups * 2 selected)


def test_create_batches_with_group(tmp_path):
    from mating_kernel.problems.mk_suite import MKSuiteConfig

    mk_config = MKSuiteConfig(cobi=True, re=True, constrained=False, num_objectives=[4])
    settings = {"setting1": [1, 2], "setting2": ["a", "b", "c"]}
    batches = Batcher.create_batches(
        suite_config=mk_config,
        settings=settings,
        reps=2,
        num_batches=None,  # Let it group by sid
        seed=123,
        exp_dir=tmp_path,
        filter_bids=None,
        skip_existing=False,
        group_by=["rep"],  # Group by rep to create 2 groups (rep 0 and rep 1)
    )
    # There should be 2 unique sids (0 and 1), so we expect 2 groups
    assert len(batches) == 2
    total_batches = sum(len(batch_group) for batch_group in batches)
    assert total_batches == 24  # Total number of batches remains the same


def test_create_batches_group_num_conflict(tmp_path):
    from mating_kernel.problems.mk_suite import MKSuiteConfig

    mk_config = MKSuiteConfig(cobi=True, re=True, constrained=False, num_objectives=[4])
    settings = {"setting1": [1, 2], "setting2": ["a", "b", "c"]}
    with pytest.raises(ValueError):
        Batcher.create_batches(
            suite_config=mk_config,
            settings=settings,
            reps=2,
            num_batches=2,  # This conflicts with group_by
            seed=123,
            exp_dir=tmp_path,
            filter_bids=None,
            skip_existing=False,
            group_by=["sid"],  # Group by sid to create 6 groups
        )
