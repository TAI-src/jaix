from jaix.runner.ask_tell.strategy import CMAConfig, CMA
from . import DummyEnv, loop
import pytest


def test_init():
    opt = CMA(CMAConfig(sigma0=0.2), [[0, 0, 0]])
    assert opt.sigma0 == 0.2
    assert opt.N_pheno == 3


@pytest.mark.parametrize("dimension,num_objectives", [(3, 1), (3, 2)])
def test_loop(dimension, num_objectives):
    opt = CMA(CMAConfig(sigma0=0.2), [[0] * dimension])

    if num_objectives > 1:
        with pytest.raises(AssertionError):
            X, Y = loop(dimension, num_objectives, opt)
    else:
        X, Y = loop(dimension, num_objectives, opt)
        assert len(X) == opt.popsize
