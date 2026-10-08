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
import pandas as pd


def test_parsing(tmp_path):
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
            str(tmp_path / "test"),
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
    assert config.out_dir == tmp_path / "test"
    assert config.settings == {
        "perf_stats_dir": ["this"],
        "kernel": [True],
        "abs": [True],
        "rel_fit": [True],
        "state": [True],
        "age": [True],
        "keep_mutated": [True],
        "no_x": [False],
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
        no_x=False,
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
    "kernel_f": {
        "kernel": True,
        "abs": False,
        "rel_fit": False,
        "state": False,
        "age": False,
        "keep_mutated": False,
        "target": "o_F_0",
        "feature_analysis": False,
        "no_x": False,
    },
    "rel_f": {
        "kernel": True,
        "abs": False,
        "rel_fit": True,
        "state": False,
        "age": True,
        "keep_mutated": False,
        "target": "o_F_0",
        "feature_analysis": False,
        "no_x": True,
    },
}


def make_batches(tmp_path, settings, num_batches=2):
    problems = [
        make_tracked(REProblem)(REProblemConfig(), inst=i) for i in range(num_batches)
    ]
    batches = [
        Batch(
            problem=problem,
            seed=123,
            pid=ProblemInfo(problem).uuid,
            sid=0,
            pinfo=ProblemInfo(problem),
            parent_dir=tmp_path,
            perf_stats_dir=str(Path(__file__).parent.parent / "data"),
            **settings,
        )
        for problem in problems
    ]
    return batches


@pytest.mark.parametrize("setting_name", list(common_settings.keys()))
def test_run_batches(tmp_path, setting_name):
    settings = common_settings[setting_name]
    batches = make_batches(tmp_path, settings)
    result_files = PredExperiment._run_batches(batches)
    assert len(result_files) == len(batches)
    for batch, files in zip(batches, result_files):
        assert len(files) == len(PredExperiment.file_paths(batch))
        for f in files:
            assert Path(f).exists()


@pytest.mark.parametrize("setting_name", list(common_settings.keys()))
def test_post_process_batches(tmp_path, setting_name):
    settings = common_settings[setting_name]
    batches = make_batches(tmp_path, settings)

    PredExperiment._run_batches(batches)
    post_files = PredExperiment._post_process_batches(batches)
    assert len(post_files) == 1
    assert len(post_files[0]) == 4
    for f in post_files[0]:
        assert Path(f).exists()
        if f.endswith(".csv"):
            data = pd.read_csv(f)  # Check if the file can be read as a CSV
            assert not data.empty  # Check if the DataFrame is not empty
            assert "cv_score_mean" in data.columns  # Check for expected column
            assert "cv_score_std" in data.columns  # Check for expected column
            assert "setting" in data.columns  # Check for expected column
            assert "pid" in data.columns  # Check for expected column
            assert "RE22" in data["pid"].values  # Check for expected pid value
