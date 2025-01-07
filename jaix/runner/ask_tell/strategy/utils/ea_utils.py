import numpy as np
from typing import List


def global_flip(parent, p=None, low=0, high=1):
    x = parent.copy()
    if p is None:
        p = 1 / len(x)
    assert p >= 0 and p <= 1
    for i in range(len(x)):
        if np.random.rand() < p:
            options = list(range(low, high + 1))
            options.remove(x[i])
            x[i] = np.random.choice(options)
    return x


def onepoint_crossover(x1, x2, k=None):
    n = len(x1)
    assert len(x2) == n
    if k is None:
        k = np.random.randint(0, n)
    else:
        assert k >= 0 and k < n
    return np.concatenate([x1[:k], x2[k:]])


def uniform_crossover(x1, x2, mask=None):
    n = len(x1)
    assert len(x2) == n
    if mask is None:
        mask = np.random.randint(0, 2, n)
    else:
        assert len(mask) == n
        assert all([m == 0 or m == 1 for m in mask])
    return np.where(mask, x1, x2)


class Individual:
    def __init__(self, x, fitness):
        self.x = x
        self.fitness = fitness


def select(population: List[Individual], mu: int):
    return sorted(population, key=lambda x: x.fitness)[:mu]
