from .rbf.test_rbf_adapter import get_config
from jaix.env.utils.problem import RBFFitConfig, RBFFit
from jaix.env.singular import (
    ECEnvironment,
    ECEnvironmentConfig,
)
from ttex.config import ConfigurableObjectFactory as COF
import pytest


def test_rbf_fit():
    rbf_adapter_config = get_config()
    config = RBFFitConfig(rbf_adapter_config, 1e-8)
    rbf = RBFFit(config, 5)
    x = [1] * rbf.dimension
    assert isinstance(rbf._eval(x)[0], float)
    assert not rbf.final_target_hit()


# integration test env
def test_with_env():
    rbf_adapter_config = get_config()
    config = RBFFitConfig(rbf_adapter_config, 1e-8)
    func = COF.create(RBFFit, config, 10)
    config = ECEnvironmentConfig(budget_multiplier=1)
    env = COF.create(ECEnvironment, config, func)

    info = env._get_info()
    env.step(env.action_space.sample())
