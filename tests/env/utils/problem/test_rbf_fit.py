from .rbf.test_rbf_adapter import get_config
from jaix.env.utils.problem.rbf_fit import RBFFitConfig, RBFFit
from jaix.env.singular.ec_env import (
    ECEnvironment,
    ECEnvironmentConfig,
)
from ttex.config import ConfigurableObjectFactory as COF


def test_rbf_fit():
    rbf_adapter_config = get_config()
    config = RBFFitConfig(rbf_adapter_config, 1e-8)
    rbf = RBFFit(config, 5)
    x = [1] * rbf.dimension
    fit, _ = rbf._eval(x)
    assert isinstance(fit[0], float)
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


def test_instance_seeding():
    rbf_adapter_config = get_config()
    rbf_adapter_config.noisy = False
    config = RBFFitConfig(rbf_adapter_config, 1e-8)
    rbf1 = RBFFit(config, 5)
    rbf2 = RBFFit(config, 5)
    rbf3 = RBFFit(config, 6)
    x = [1] * rbf1.dimension
    assert rbf1._eval(x) == rbf2._eval(x)
    assert rbf1.dimension != rbf3.dimension or rbf1._eval(x) != rbf3._eval(x)
