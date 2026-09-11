import numpy as np


def angle_between(
    vec1: np.ndarray, vec2: np.ndarray, axis: int = -1, min_norm: float = 1e-12
) -> np.ndarray | np.floating:
    """
    Compute the angle between two vectors in radians.
    """
    norm1 = np.linalg.norm(vec1, axis=axis)
    norm2 = np.linalg.norm(vec2, axis=axis)

    valid = (norm1 > min_norm) & (norm2 > min_norm)

    dot_product = np.sum(vec1 * vec2, axis=axis)
    cos_angle = np.divide(
        dot_product,
        norm1 * norm2,
        out=np.full_like(dot_product, np.nan, dtype=float),
        where=valid,
    )

    return np.arccos(np.clip(cos_angle, -1.0, 1.0))
