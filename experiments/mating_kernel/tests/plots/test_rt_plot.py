from pathlib import Path

import pandas as pd
import pytest
from cocoviz import ProblemDescription, Result, ResultSet
from cocoviz.exceptions import NoSuchIndicatorException
from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from mating_kernel.plots.rt_plot import (
    get_description,
    get_result_set,
    plot_rt,
    plot_rt_cocoviz,
    results_from_df,
)
from mating_kernel.problems.cobi_configs import get_config
from mating_kernel.problems.problem_info import ProblemInfo

test_folder = Path(__file__).parent.parent / "data" / "cocoviz"


cobi_problem = CobiProblem(get_config(0), inst=0)
re_problem = REProblem(REProblemConfig(), 0)


test_data_regular = [
    {
        "a_coverage": 0.5,
        "a_score": 0.6,
        "a_fevals": 200,
        "seed": 1,
    },
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


@pytest.mark.parametrize("problem", [cobi_problem, re_problem])
def test_get_description(problem):
    problem_info = ProblemInfo(problem)
    desc = get_description(problem_info)  # Ensure it runs without errors
    assert desc.name == problem_info.uuid
    assert desc.number_of_variables == problem_info.num_variables
    assert isinstance(desc, ProblemDescription)


def test_from_df():
    df = pd.DataFrame(test_data_regular)
    res_set = results_from_df(
        df,
        "test_exp",
        ProblemDescription("test_problem", 0, 10, 2),
        ["a_coverage", "a_score"],
    )
    assert (
        len(res_set) == df["seed"].nunique()
    )  # Check that the number of results matches the number of unique seeds
    assert len(res_set.algorithms) == 1
    for result in res_set:
        assert isinstance(result, Result)
        assert hasattr(result, "pq_file_name")
        assert (
            f"_{result.seed}.parquet" in result.pq_file_name
        )  # Check that the seed is in the pq_file_name


@pytest.mark.parametrize("with_path", [True, False])
def test_get_result_set_no_df(tmp_path, with_path):
    params = {
        "source_df": None,
        "exp_name": "test_exp",
        "pdesc": ProblemDescription("test_problem", 0, 10, 2),
        "indicators": ["a_coverage", "a_score"],
        "pq_path": tmp_path if with_path else None,
    }
    if not with_path:
        with pytest.raises(ValueError, match="source_df is None and pq_path is None"):
            get_result_set(**params)
    else:
        result_set = get_result_set(**params)
        assert isinstance(result_set, ResultSet)
        assert (
            len(result_set) == 0
        )  # Since there are no existing files, the result set should be empty


@pytest.mark.parametrize("with_path", [True, False])
def test_get_result_set_with_df(tmp_path, with_path):
    df = pd.DataFrame(test_data_regular)
    params = {
        "source_df": df,
        "exp_name": "test_exp",
        "pdesc": ProblemDescription("test_problem", 0, 10, 2),
        "indicators": ["a_coverage", "a_score"],
        "pq_path": tmp_path if with_path else None,
    }
    result_set = get_result_set(**params)
    assert isinstance(result_set, ResultSet)
    assert len(result_set) == df["seed"].nunique()
    if with_path:
        # Check that the parquet files are created
        for result in result_set:
            assert (tmp_path / result.pq_file_name).exists()

        # also check that we can load it back from the directory
        params["source_df"] = None
        result_set_loaded = get_result_set(**params)
        assert isinstance(result_set_loaded, ResultSet)
        assert len(result_set_loaded) == len(result_set)


def test_plot_rt(tmp_path):
    # load data directry from test data folder
    res_set = ResultSet.from_directory(test_folder)
    indicator_specs = [("a_score", True), ("a_coverage", True)]
    file_names = plot_rt(
        res_set, indicator_specs=indicator_specs, output_dir=tmp_path, file_prefix="a"
    )
    assert len(file_names) == len(indicator_specs)
    for plot_path in file_names:
        # check that the file exists
        assert Path(plot_path).exists()


@pytest.mark.parametrize("with_path", [True, False])
def test_plot_rt_cocoviz_df(tmp_path, with_path):
    if not with_path:
        src = pd.DataFrame(test_data_regular)
    else:
        src = test_folder
    file_names = plot_rt_cocoviz(
        source=src,
        exp_name="test",
        pinfo=ProblemInfo(cobi_problem),
        output_dir=tmp_path,
        indicator_specs=[("a_score", True), ("a_coverage", True)],
        file_prefix="",
    )
    assert len(file_names) == 2


def test_default_indicators(tmp_path):
    with pytest.raises(NoSuchIndicatorException):
        file_names = plot_rt_cocoviz(
            source=test_folder,
            exp_name="test",
            pinfo=ProblemInfo(cobi_problem),
            output_dir=tmp_path,
            indicator_specs=None,
            file_prefix="",
        )
