from pathlib import Path

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


def test_parsing():
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
            "results",
            "--nth",
            "0",
            "1",
            "--selector",
            "random",
            "default",
            "--n_gen",
            "1000",
        ]
    )
    assert config.suite_config.cobi is True
    assert config.reps == 5
    assert config.num_batches == 1
    assert config.seed == 123
    assert config.out_dir == Path("results")
    assert config.settings == {
        "alg_name": ["NSGA2"],
        "selector": ["random", "default"],
        "n_gen": [1000],
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


def test_get_recorded_alg():
    algorithm = PerfExperiment.get_recorded_alg(
        algorithm_name="NSGA2",
        selector=None,
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


def test_run_instrumented_pymoo():
    tracked_REProblem = make_tracked(REProblem)
    problem = tracked_REProblem(REProblemConfig(), inst=0)
    result, record_stats = PerfExperiment.run_instrumented_pymoo(
        problem=problem,
        algorithm_name="NSGA2",
        selector=None,
        n_gen=5,
        algorithm_params={"pop_size": 10},
        seed=123,
    )
    assert result is not None
    assert isinstance(record_stats, pd.DataFrame)


def test_run_batch(tmp_path):
    problem = make_tracked(REProblem)(REProblemConfig(), inst=0)
    pinfo = ProblemInfo(problem)
    batch = Batch(
        name="test_batch",
        problem=problem,
        alg_name="NSGA2",
        selector=None,
        n_gen=5,
        seed=123,
        pid=pinfo.uuid,
        sid=0,
        pinfo=pinfo,
        parent_dir=tmp_path,
    )
    batch_cpy = batch.model_copy()
    result_files = PerfExperiment._run_batch(batch)
    assert len(result_files) == 3
    result_file, record_stats_file, batch_file = result_files
    assert Path(result_file).exists()
    assert Path(record_stats_file).exists()
    assert Path(batch_file).exists()
    # check that the batch and problem are not modified
    assert batch == batch_cpy
    assert batch.problem == batch_cpy.problem
