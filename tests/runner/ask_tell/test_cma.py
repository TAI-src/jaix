from jaix.runner.ask_tell.strategy.cma import CMAConfig, CMA
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


@pytest.mark.parametrize("warm_start_best", [True, False])
def test_warm_start(warm_start_best):
    dimension = 3
    env = DummyEnv(dimension=dimension, num_objectives=1)
    opt = CMA(CMAConfig(sigma0=0.2, warm_start_best=warm_start_best), [[0] * dimension])
    for _ in range(10):
        X, Y = loop(dimension, 1, opt, env)
    assert all(opt.x0 == [0] * dimension)
    # Note: must be same dimension
    env.reset()
    xlast = X[-1]
    xfav = opt.result.xfavorite
    opt.warm_start(xlast, env)
    for _ in range(10):
        X, Y = loop(dimension, 1, opt, env)
    if warm_start_best:
        assert all(opt.x0 == xfav)
    else:
        assert all(opt.x0 == xlast)
