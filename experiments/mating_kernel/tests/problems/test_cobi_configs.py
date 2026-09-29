from mating_kernel.problems.cobi_configs import get_configs

from jaix.env.utils.problem.cobi_problem import CobiProblemConfig


def test_get_configs():
    configs = get_configs()
    assert isinstance(configs, list)
    assert len(configs) == 7
    for config in configs:
        assert isinstance(config, CobiProblemConfig)
