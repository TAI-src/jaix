from jaix.env.utils.problem.rbf.rbf import RBFKernel, RBF
import numpy as np


def test_kernel():
    kernel = RBFKernel.GAUSSIAN
    assert kernel(0, 5) == 1
    assert kernel(1, 1) == np.exp(-1)


def test_rbf_single():
    rbf = RBF([0], [1], [1], RBFKernel.GAUSSIAN)
    assert rbf.eval(0) == 1
    assert rbf.eval(1) == np.exp(-1)


def test_rbf_double():
    rbf = RBF([0, 1], [1, 1], [0.2, 0.7], RBFKernel.GAUSSIAN)
    assert rbf.eval(0) == 0.2 * 1 + 0.7 * np.exp(-1)
    assert rbf.eval(1) == 0.7 * 1 + 0.2 * np.exp(-1)

    rbf = RBF([0, 1], [0, 1], [0.2, 0.7], RBFKernel.GAUSSIAN)
    assert rbf.eval(0) == 0.2 + 0.7 * np.exp(-1)
    assert rbf.eval(1) == 0.2 + 0.7 * 1
