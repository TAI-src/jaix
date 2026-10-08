from pathlib import Path
import pytest
import pandas as pd
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)

from mating_kernel.exp.perf_exp import PerfExperiment
from mating_kernel.exp.utils.batch import Batch
from mating_kernel.problems.mo_tracking import make_tracked
from mating_kernel.problems.problem_info import ProblemInfo
from mating_kernel.pymoo.parser.reproduction_parser import ReproductionParser
from mating_kernel.pymoo.recording_callback import recursive_getattr


def test_parsing(tmp_path):
    # Test parser together with the parse base experiment
    config, filter_bids = PerfExperiment.parse_args(
        [
            "--cobi",
            "--reps",
            "5",
            "--num_batches",
            "1",
            "--seed",
            "123",
            "--out_dir",
            str(tmp_path / "results"),
            "--nth",
            "0",
            "1",
            "--selector",
            "random",
            "default",
            "--n_gen",
            "1000",
            "--num_candidates",
            "3",
            "1",
        ]
    )
    assert config.suite_config.cobi is True
    assert config.reps == 5
    assert config.num_batches == 1
    assert config.seed == 123
    assert config.out_dir == tmp_path / "results"
    assert config.settings == {
        "alg_name": ["NSGA2"],
        "selector": ["random", "default"],
        "n_gen": [1000],
        "num_candidates": [3, 1],
        "candidate_pressure": [2],
    }
    assert filter_bids == [0, 1]


def test_get_alg_class():
    alg_class = PerfExperiment.get_alg_class("NSGA2")
    from pymoo.algorithms.moo.nsga2 import NSGA2

    assert alg_class == NSGA2

    try:
        PerfExperiment.get_alg_class("UnsupportedAlg")
    except ValueError as e:
        assert str(e) == "Unsupported algorithm: UnsupportedAlg"
    else:
        raise AssertionError("Expected get_alg_class to reject unsupported algorithm")


@pytest.mark.parametrize("selector", ["default", "random"])
def test_get_recorded_alg(selector):
    algorithm = PerfExperiment.get_recorded_alg(
        algorithm_name="NSGA2",
        selector=selector,
        selector_params={"num_candidates": 3, "candidate_pressure": 2},
        algorithm_params={"pop_size": 5},
        record_args=ReproductionParser.record_args,
        record_attributes=ReproductionParser.record_attributes,
    )
    from pymoo.algorithms.moo.nsga2 import NSGA2

    assert isinstance(algorithm, NSGA2)
    assert algorithm.pop_size == 5
    record_ret = ReproductionParser.record_retrieval
    for attr in record_ret:
        attr = recursive_getattr(algorithm, attr)
        # if the value has retrieve_records method, this means that it has been instrumented correctly
        # and we will be able to retrieve the values later
        assert attr is not None and hasattr(attr, "retrieve_records")


@pytest.mark.parametrize("selector", ["default", "random"])
def test_run_instrumented_pymoo(selector):
    tracked_REProblem = make_tracked(REProblem)
    problem = tracked_REProblem(REProblemConfig(), inst=0)
    result, record_stats = PerfExperiment.run_instrumented_pymoo(
        problem=problem,
        algorithm_name="NSGA2",
        selector=selector,
        selector_params={"num_candidates": 3, "candidate_pressure": 2},
        n_gen=5,
        algorithm_params={"pop_size": 10},
        seed=123,
    )
    assert result is not None
    assert isinstance(record_stats, pd.DataFrame)


@pytest.mark.parametrize("selector", ["default", "random"])
def test_run_batch(tmp_path, selector):
    problem = make_tracked(REProblem)(REProblemConfig(), inst=0)
    pinfo = ProblemInfo(problem)
    batch = Batch(
        name="test_batch",
        problem=problem,
        alg_name="NSGA2",
        selector=selector,
        n_gen=5,
        seed=123,
        pid=pinfo.uuid,
        sid=0,
        pinfo=pinfo,
        parent_dir=tmp_path,
        num_candidates=3,
        candidate_pressure=2,
    )
    batch_cpy = batch.model_copy()
    result_files = PerfExperiment._run_batch(batch)
    assert len(result_files) == 3
    result_file, record_stats_file, archive_file = result_files
    assert Path(result_file).exists()
    assert Path(record_stats_file).exists()
    assert Path(archive_file).exists()
    # check that the batch and problem are not modified
    assert batch == batch_cpy
    assert batch.problem == batch_cpy.problem


@pytest.mark.parametrize("selector", ["default", "random"])
def test_post_process_batches(tmp_path, selector):
    problem = make_tracked(REProblem)(REProblemConfig(), inst=0)
    pinfo = ProblemInfo(problem)
    batches = []
    for i in range(3):
        batch = Batch(
            problem=problem,
            alg_name="NSGA2",
            selector=selector,
            n_gen=3,
            seed=123 + i,
            pid=pinfo.uuid,
            sid=0,
            rep=i,
            pinfo=pinfo,
            parent_dir=tmp_path,
            num_candidates=3,
            candidate_pressure=2,
        )
        PerfExperiment._run_batch(batch)
        batches.append(batch)
    merged_files = PerfExperiment._post_process_batches(batches)
    assert len(merged_files) == 1
    merged_file = merged_files[0][0]
    assert Path(merged_file).exists()
    merged_df = pd.read_csv(merged_file)
    assert (
        len(merged_df) >= 3 * 1 * 100
    )  # 3 batches * 2 generation * 100 (at least) offspring per generation
    assert "seed" in merged_df.columns
    unique_seeds = merged_df["seed"].unique()
    assert set(unique_seeds) == {123, 124, 125}  # seeds used in the batches
