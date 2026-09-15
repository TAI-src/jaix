import numpy as np
from pymoo.algorithms.moo.nsga3 import (
    HyperplaneNormalization,
    ReferenceDirectionSurvival,
)
from pymoo.core.survival import Survival


class StaticHyperplaneNormalization(HyperplaneNormalization):
    def __init__(self, ideal_point: np.ndarray, nadir_point: np.ndarray):
        super().__init__(n_dim=len(ideal_point))
        self.ideal_point = ideal_point
        self.nadir_point = nadir_point

    def update(self, F, nds=None):
        pass  # Do not update the normalization, keep it static


class StaticReferenceDirectionSurvival(ReferenceDirectionSurvival):
    def __init__(
        self, ref_dirs: np.ndarray, ideal_point: np.ndarray, nadir_point: np.ndarray
    ):
        Survival.__init__(self, filter_infeasible=True)
        self.ref_dirs = ref_dirs
        self.opt = None
        self.norm = StaticHyperplaneNormalization(ideal_point, nadir_point)
