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


def test_seed_batches(tmp_path):
    batch = make_batch(tmp_path)

    result = Batcher.seed_batches(
        [batch],
        reps=3,
        seed=123,
    )

    assert len(result) == 3
    assert [b.rep for b in result] == [0, 1, 2]
    assert all(b.seed is not None for b in result)
    # simulated run experiments by creating files in outdir
    for b in result:
        sim_run(b)
    # remove out_dir for batch2 and see that it is picked up again
    seed1 = result[1].seed
    shutil.rmtree(result[1].out_dir)
    result1 = Batcher.seed_batches(
        [batch],
        reps=3,
        seed=123,
    )
    assert len(result1) == 1
    assert result1[0].rep == 1
    assert result1[0].seed == seed1

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
