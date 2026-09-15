import numpy as np
from pymoo.algorithms.moo.nsga3 import NSGA3, ReferenceDirectionSurvival
from pymoo.core.population import Population

from utils_nsga3_norm import (
    StaticHyperplaneNormalization,
    StaticReferenceDirectionSurvival,
)


def test_static_hyperplane_normalization():
    static_norm = StaticHyperplaneNormalization(
        ideal_point=np.array([0.0, 0.0]), nadir_point=np.array([1.0, 1.0])
    )

    static_norm.update(F=np.array([[0.5, 0.5], [0.2, 0.8]]))
    assert np.array_equal(static_norm.ideal_point, np.array([0.0, 0.0]))
    assert np.array_equal(static_norm.nadir_point, np.array([1.0, 1.0]))


def test_static_reference_direction_survival():
    static_survival = StaticReferenceDirectionSurvival(
        ref_dirs=np.array([[1.0, 0.0], [0.0, 1.0]]),
        ideal_point=np.array([0.0, 0.0]),
        nadir_point=np.array([1.0, 1.0]),
    )
    pop = Population.new("X", np.array([[0.5, 0.5], [0.2, 0.8]]))
    pop.set("F", np.array([[0.5, 0.5], [0.2, 0.8]]))
    survived = static_survival._do(None, pop, 1)
    assert len(survived) == 1
    assert np.array_equal(static_survival.norm.ideal_point, np.array([0.0, 0.0]))
    assert np.array_equal(static_survival.norm.nadir_point, np.array([1.0, 1.0]))

    # test normal survival with updated normalization
    survival = ReferenceDirectionSurvival(ref_dirs=np.array([[1.0, 0.0], [0.0, 1.0]]))
    survived = survival._do(None, pop, 1)
    assert not np.array_equal(survival.norm.ideal_point, np.array([0.0, 0.0]))
    assert not np.array_equal(survival.norm.nadir_point, np.array([1.0, 1.0]))


def test_within_nsga3():
    survival = StaticReferenceDirectionSurvival(
        ref_dirs=np.array([[1.0, 0.0], [0.0, 1.0]]),
        ideal_point=np.array([0.0, 0.0]),
        nadir_point=np.array([1.0, 1.0]),
    )
    algorithm = NSGA3(
        pop_size=2, ref_dirs=np.array([[1.0, 0.0], [0.0, 1.0]]), survival=survival
    )
    assert algorithm.survival.norm.ideal_point is not None
