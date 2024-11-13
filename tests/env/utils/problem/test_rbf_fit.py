from .rbf.test_rbf_adapter import get_config
from jaix.env.utils.problem import RBFFitConfig, RBFFit

def test_rbf_fit():
    rbf_adapter_config = get_config()
    config = RBFFitConfig(rbf_adapter_config, 1e-8)
    rbf = RBFFit(config, 5)
    x = [0]*10
    print(rbf._eval(x))
    # TODO: proper test

