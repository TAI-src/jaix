import pytest
from jaix.env.utils.archive.mo_archive import MOArchiveConfig
from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)

from nsga3_experiment import NSGA3ExperimentConfig


def test_create_mo_archive_config():
    problem = REProblem(REProblemConfig(), inst=2)
    config = NSGA3ExperimentConfig.create_mo_archive_config(problem)
    assert config.max_size is not None
    assert isinstance(config.num_refpoints, int)
    assert isinstance(config, MOArchiveConfig)


def test_re_problem_list():
    problems = NSGA3ExperimentConfig.re_problem_list()
    assert isinstance(problems, list)
    assert len(problems) > 0
    assert all(isinstance(p, REProblem) for p in problems)


def test_cobi_problem_list():
    problems = NSGA3ExperimentConfig.cobi_problem_list()
    assert isinstance(problems, list)
    assert len(problems) > 0
    assert all(isinstance(p, CobiProblem) for p in problems)


@pytest.mark.parametrize("mode", ["cobi", "reproblem", ""])
def test_generate_problem_list(mode):
    problems = NSGA3ExperimentConfig.generate_problem_list(mode)
    assert isinstance(problems, list)
    assert len(problems) > 0
    if mode == "cobi":
        assert all(isinstance(p, CobiProblem) for p in problems)
    elif mode == "reproblem":
        assert all(isinstance(p, REProblem) for p in problems)
    else:
        assert all(isinstance(p, (CobiProblem, REProblem)) for p in problems)


@pytest.mark.parametrize("kwargs", [{"max_size": 10}, None, {"num_refpoints": 50}])
def test_config_init(kwargs):
    config = NSGA3ExperimentConfig(
        num_independent_runs=1,
        num_generations=10,
        mo_archive_kwargs=kwargs,
    )
    assert config.num_prefill_samples == -1
    problem = REProblem(REProblemConfig(), inst=2)
    config.update_defaults(problem)
    assert isinstance(config.mo_archive_config, MOArchiveConfig)
    assert config.mo_archive_config.max_size is not None
    if kwargs is not None and "max_size" in kwargs:
        assert config.mo_archive_config.max_size == kwargs["max_size"]
    if kwargs is not None and "num_refpoints" in kwargs:
        assert config.mo_archive_config.num_refpoints == kwargs["num_refpoints"]
    assert config.num_prefill_samples == config.mo_archive_config.num_refpoints
    assert config.num_offspring == config.mo_archive_config.num_refpoints
