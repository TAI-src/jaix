from utils_space import angle_between
import numpy as np


def test_angle_between():
    # Test with 1D vectors
    vec1 = np.array([1, 0])
    vec2 = np.array([0, 1])
    angle = angle_between(vec1, vec2)
    assert np.isclose(angle, np.pi / 2)

    # Test with 2D vectors
    vec1 = np.array([[1, 0], [0, 1]])
    vec2 = np.array([[0, 1], [1, 0]])
    angles = angle_between(vec1, vec2)
    assert np.allclose(angles, [np.pi / 2, np.pi / 2])

    # Test with zero vector
    vec1 = np.array([0, 0])
    vec2 = np.array([1, 0])
    angle = angle_between(vec1, vec2)
    assert np.isnan(angle)

    # Test with min_norm threshold
    vec1 = np.array([1e-13, 0])
    vec2 = np.array([1, 0])
    angle = angle_between(vec1, vec2)
    assert np.isnan(angle)
