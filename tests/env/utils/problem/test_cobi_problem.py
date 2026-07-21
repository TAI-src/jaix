from . import CobiProblem, CobiProblemConfig
import numpy as np
import pytest


@pytest.fixture(scope="session", autouse=True)
def skip_remaining_tests():
    try:
        import cobi  # noqa: F401

        assert CobiProblem is not None
    except ImportError:
        assert CobiProblem is None
        pytest.skip(
            "Skipping CobiProblem tests. If this is unexpected, check that the jaix_cobi docker image is used."
        )


def test_cobi_problem():
    config = CobiProblemConfig(n_var=3)
    problem = CobiProblem(config, 1)
    f, _ = problem([0, 0, 0])
    assert len(f) == problem.num_objectives
    assert problem.lower_bounds == [config.domain[0]] * config.n_var
    assert problem.upper_bounds == [config.domain[1]] * config.n_var


def test_cobi_feasible():
    config = CobiProblemConfig(
        n_var=3, n_constraints={"Linear": 0, "Quadratic": 0, "Multi": 0}, domain=(-4, 4)
    )
    problem = CobiProblem(config, 1)
    x_feasible = [0, 0, 0]
    f, _ = problem(x_feasible)
    assert not any(np.isnan(f))  # Check that the feasible point returns valid values

    x_infeasible = [5, 5, 5]  # Outside the domain, should be infeasible
    f, _ = problem(x_infeasible)
    assert all(np.isnan(f))  # Check that the infeasible point returns NaN values
