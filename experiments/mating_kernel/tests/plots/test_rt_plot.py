from pathlib import Path

from cocoviz import ResultSet
from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from mating_kernel.plots.rt_plot import (
    plot_rt,
)
from mating_kernel.problems.cobi_configs import get_config

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
