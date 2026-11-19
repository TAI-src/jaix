from mating_kernel.problems.problem_info import ProblemInfo
from mating_kernel.exp.utils.batch import Batch
from mating_kernel.plots.cv_result import Result, ResultSet
from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
import pytest
from mating_kernel.problems.cobi_configs import get_config
from cocoviz import ProblemDescription
import pandas as pd

cobi_problem = CobiProblem(get_config(0), inst=0)
re_problem = REProblem(REProblemConfig(), 0)


@pytest.mark.parametrize("problem", [cobi_problem, re_problem])
def test_get_description(problem):
    problem_info = ProblemInfo.from_problem(problem)
    desc = Result.get_description(problem_info)  # Ensure it runs without errors
    assert desc.name == problem_info.uuid
    assert desc.number_of_variables == problem_info.num_variables
    assert isinstance(desc, ProblemDescription)


def create_test_batch(problem):
    problem_info = ProblemInfo.from_problem(problem)
    return Batch(
        pid="test_pid",
        sid=1,
        pinfo=problem_info,
        parent_dir=".",
        seed=42,
    )


test_data = [
    {
        "a_coverage": 0.5,
        "a_score": 0.6,
        "a_fevals": 200,
        "seed": 1,
    },
    {
        "a_coverage": 0.6,
        "a_score": 0.7,
        "a_fevals": 300,
        "seed": 1,
    },
    {
        "a_coverage": 0.7,
        "a_score": 0.8,
        "a_fevals": 200,
        "seed": 2,
    },
]


def test_init_result():
    batch = create_test_batch(cobi_problem)
    result = Result(data=pd.DataFrame(test_data), batch=batch)
    assert result.algorithm == "_".join(f"{k}={v}" for k, v in batch.settings.items())
    assert result.problem.name == batch.pinfo.uuid


def test_write_read(tmp_path):
    batch = create_test_batch(cobi_problem)
    result = Result(data=pd.DataFrame(test_data), batch=batch)
    written_path = result.write(tmp_path)
    expected_file = tmp_path / f"{batch.run_id}.parquet"
    assert written_path == expected_file
    assert expected_file.exists()

    # read back the result and check if it matches
    res_read = Result.from_parquet(expected_file)
    assert res_read.algorithm == result.algorithm
    assert res_read.problem.name == result.problem.name
    assert res_read.batch.pid == result.batch.pid
    assert res_read.batch.sid == result.batch.sid
    assert res_read.batch.seed == result.batch.seed


def test_result_set(tmp_path):
    batch1 = create_test_batch(cobi_problem)
    result1 = Result(data=pd.DataFrame(test_data), batch=batch1)
    result1.write(tmp_path)
    batch2 = create_test_batch(re_problem)
    result2 = Result(data=pd.DataFrame(test_data), batch=batch2)
    # writing to subdir to test that ResultSet can handle different directories
    result2.write(tmp_path / "subdir")
    res_set = ResultSet.from_directory(tmp_path)
    assert len(res_set) == 2
    for res in res_set:
        assert isinstance(res, Result)
        assert res.batch is not None
        assert isinstance(res.batch, Batch)
        assert hasattr(res, "pq_file_name")
