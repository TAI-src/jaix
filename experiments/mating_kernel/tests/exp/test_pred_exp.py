from mating_kernel.exp.pred_exp import PredExperiment
from pathlib import Path
from mating_kernel.exp.utils.batch import Batch

from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)

from mating_kernel.problems.mo_tracking import make_tracked
from mating_kernel.problems.problem_info import ProblemInfo
import pytest


def test_parsing():
    config, filter_bids = PredExperiment.parse_args(
        [
            "--cobi",
            "--reps",
            "5",
            "--num_batches",
            "1",
            "--seed",
            "123",
            "--out_dir",
            "test",
            "--nth",
            "0",
            "1",
            "--perf_stats_dir",
            "this",
            "--kernel",
            "--abs",
            "--rel_fit",
            "--state",
            "--age",
            "--keep_mutated",
            "--target",
            "survived",
            "--feature_analysis",
        ]
    )
    assert config.suite_config.cobi is True
    assert config.reps == 5
    assert config.num_batches == 1
    assert config.seed == 123
    assert config.out_dir == Path("test")
    assert config.settings == {
        "perf_stats_dir": ["this"],
        "kernel": [True],
        "abs": [True],
        "rel_fit": [True],
        "state": [True],
        "age": [True],
        "keep_mutated": [True],
        "target": ["survived"],
        "feature_analysis": [True],
    }
    assert filter_bids == [0, 1]


def test_get_batch(tmp_path):
    problem = make_tracked(REProblem)(REProblemConfig(), inst=0)
    pinfo = ProblemInfo(problem)
    batch = Batch(
        problem=problem,
        n_gen=5,
        seed=123,
        pid=pinfo.uuid,
        sid=0,
        pinfo=pinfo,
        parent_dir=tmp_path,
        perf_stats_dir=str(Path(__file__).parent.parent / "data"),
        kernel=True,
        abs=True,
        rel_fit=True,
        state=True,
        age=True,
        keep_mutated=True,
        target="survived",
        feature_analysis=False,
    )
    batch_cpy = batch.model_copy()
    result_files = PredExperiment._run_batch(batch)
    assert len(result_files) == len(PredExperiment.file_paths(batch))
    for f in result_files:
        assert Path(f).exists()
    assert batch == batch_cpy
    assert batch.problem == batch_cpy.problem


common_settings = {
    "f_pred": {
        "kernel": True,
        "abs": False,
        "rel_fit": False,
        "state": False,
        "age": False,
        "keep_mutated": False,
        "target": "o_F_0",
        "feature_analysis": False,
    }
}


@pytest.mark.parametrize("setting_name", list(common_settings.keys()))
def test_run_batches(tmp_path, setting_name):
    settings = common_settings[setting_name]
    problem = make_tracked(REProblem)(REProblemConfig(), inst=0)
    pinfo = ProblemInfo(problem)
    batches = [
        Batch(
            problem=problem,
            n_gen=5,
            seed=123,
            pid=pinfo.uuid,
            sid=i,
            pinfo=pinfo,
            parent_dir=tmp_path,
            perf_stats_dir=str(Path(__file__).parent.parent / "data"),
            **settings
        )
        for i in range(2)
    ]
    result_files = PredExperiment._run_batches(batches)
    assert len(result_files) == len(batches)
    for batch, files in zip(batches, result_files):
        assert len(files) == len(PredExperiment.file_paths(batch))
        for f in files:
            assert Path(f).exists()
