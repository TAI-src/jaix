from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from jaix.env.singular.ec_env import ECEnvironment, ECEnvironmentConfig
from ttex.config import ConfigurableObjectFactory as COF
import numpy as np
import pytest


@pytest.mark.parametrize("inst", list(range(23)))
def test_reproblem_init(inst):
    problem = REProblem(REProblemConfig(), inst)
    if inst >= 16:
        assert problem.constrained
    else:
        assert not problem.constrained

    assert problem.lower_bounds[0] != -np.inf  # This would be the default setting
    assert problem.evaluations == 0
    assert problem.dimension > 0
    assert problem.num_objectives > 1

    # Look at the ideal and nadir points
    assert problem.ideal_point is not None
    assert problem.ideal_point.shape[0] == problem.num_objectives
    assert problem.nadir_point is not None
    assert problem.nadir_point.shape[0] == problem.num_objectives


def test_invalid_instance():
    with pytest.raises(ValueError):
        COF.create(REProblem, REProblemConfig(), 100)
    with pytest.raises(ValueError):
        COF.create(REProblem, REProblemConfig(), -1)


@pytest.mark.parametrize("inst", list(range(23)))
def test_eval(inst):
    problem = COF.create(REProblem, REProblemConfig(), inst)
    x = np.random.uniform(problem.lower_bounds, problem.upper_bounds)
    y_raw, y_noise = problem(x)
    assert len(y_raw) == problem.num_objectives
    # check that the values are either all np.nan or numeric
    assert all(np.isnan(y) for y in y_raw) or all(np.isfinite(y) for y in y_raw)
    print(f"Instance {inst}: x={x}, y_raw={y_raw}, y_noise={y_noise}")
    assert problem.evaluations == 1


def test_env_integration():
    config = REProblemConfig()
    func = COF.create(REProblem, config, 1)
    env_config = ECEnvironmentConfig(budget_multiplier=1)
    env = COF.create(ECEnvironment, env_config, func, 0, 1)

    for _ in range(3):  # Test a few steps
        random_point = env.action_space.sample()
        obs, reward, _, truncated, info = env.step(random_point)
        assert len(obs) == func.num_objectives
        assert isinstance(info, dict)
        assert reward is None
    assert truncated
