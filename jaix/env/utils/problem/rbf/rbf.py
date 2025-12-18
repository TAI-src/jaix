from typing import List
from enum import Enum

# Based on https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RBFInterpolator.html#scipy.interpolate.RBFInterpolator
import numpy as np


class RBFKernel(Enum):
    """
    ‘linear’ : -r
    ‘thin_plate_spline’ : r**2 * log(r)
    ‘cubic’ : r**3
    ‘quintic’ : -r**5
    ‘multiquadric’ : -sqrt(1 + r**2)
    ‘inverse_multiquadric’ : 1/sqrt(1 + r**2)
    ‘inverse_quadratic’ : 1/(1 + r**2)
    ‘gaussian’ : exp(-r**2)
    """

    # LINEAR = lambda r, eps: -r
    # THIN_PLATE_SPLINE = lambda r, eps: r**2 * np.log(r)
    GAUSSIAN = lambda r, eps: np.exp(-((r * eps) ** 2))


class RBF:
    def __init__(
        self,
        c: List[float],
        eps: List[float],
        w: List[float],
        kernel: RBFKernel = RBFKernel.GAUSSIAN,
    ):
        """
        Radial Basis Function
        :param c: centers
        :param eps: shape parameters
        :param w: weights
        :param kernel: kernel function
        """
        self.c = c
        self.eps = eps
        self.w = w
        self.kernel = kernel
        assert len(c) == len(eps)
        assert len(eps) == len(w)
        # TODO: constraints on epsilon?
        # TODO: constraints on weights?

    def eval(self, x):
        vals = [
            w * self.kernel(x - c, eps) for (c, eps, w) in zip(self.c, self.eps, self.w)
        ]
        return sum(vals)
