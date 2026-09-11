from config_nsga3x import NSGA3ExperimentConfig
from jaix.env.utils.archive.mo_archive import MOArchiveConfig
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)


def test_create_mo_archive_config():
    problem = REProblem(REProblemConfig(), inst=2)
    config = NSGA3ExperimentConfig.create_mo_archive_config(problem)
    assert config.max_size is not None
    assert isinstance(config.num_refpoints, int)
    assert isinstance(config, MOArchiveConfig)
