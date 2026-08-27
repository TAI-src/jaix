from cobi_config_generator import get_configs
from jaix.env.utils.problem.cobi_problem import CobiProblemConfig


def test_get_configs():
    configs = get_configs()
    assert isinstance(configs, list)
    assert len(configs) > 0
    for config in configs:
        assert isinstance(config, CobiProblemConfig)
