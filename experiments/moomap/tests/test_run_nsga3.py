import pandas as pd
import pytest
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)

from config_run_nsga3 import get_batches
from run_nsga3 import run, run_algorithm


@pytest.mark.parametrize("static_ref", [True, False])
def test_run_algorithm(tmp_path, static_ref):
    problem = REProblem(REProblemConfig(), inst=0)
    file_path = run_algorithm(
        n_gen=2, seed=42, out_dir=tmp_path, problem=problem, static_ref=static_ref
    )
    assert file_path.endswith(".csv")
    assert tmp_path.joinpath(file_path).exists()
    assert "s42" in file_path
    assert str(problem) in file_path
    if static_ref:
        assert "fixed" in file_path  # Because static_ref=True
    else:
        assert "fixed" not in file_path  # Because static_ref=False

    df = pd.read_csv(file_path)
    assert "generation" in df.columns
    assert len(df) == 2  # Should have 2 generations
    assert "size" in df.columns


def test_run_function(tmp_path):
    args = {
        "n_gen": 2,
        "seed": None,
        "static_ref": True,
        "problem_ids": [0, 7],
        "out_dir": tmp_path,
        "n_runs": 2,
        "batch_ids": [0, 1],
        "verbose": False,
    }
    namespace = type("Args", (), args)  # Create a simple namespace object
    files = run(namespace)
    expected_batches, _ = get_batches(
        batch_id=args["batch_ids"],
        n_runs=args["n_runs"],
        problem_ids=args["problem_ids"],
        seed=args["seed"],
        static_ref=args["static_ref"],
    )
    assert len(files) == len(expected_batches)
    for batch, file in zip(expected_batches, files):
        assert str(batch["problem"]) in file
        assert "fixed" in file  # Because static_ref=True
