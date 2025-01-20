import numpy as np
from typing import List
from uuid import uuid4
import math


def global_flip(parent, p=None, low=0, high=1):
    x = parent.copy()
    if p is None:
        p = 1 / len(x)
    if np.issubdtype(type(low), int):
        low = [low] * len(x)
    if np.issubdtype(type(high), int):
        high = [high] * len(x)
    assert p >= 0 and p <= 1
    for i in range(len(x)):
        if np.random.rand() < p:
            options = list(range(low[i], high[i] + 1))
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
    def __init__(self, x, fitness, generation: int):
        self.x = x
        self.fitness = fitness
        self.generation = generation
        self.id = uuid4()

    def __repr__(self):
        return f"Individual(x={self.x}, fitness={self.fitness}, generation={self.generation})"


def select(population: List[Individual], mu: int, reverse=False):
    # TODO: test new sorting
    return sorted(
        population, key=lambda x: (x.fitness, -x.generation), reverse=reverse
    )[:mu]


def ddl_update(old_pop, new_pop, mutation_opts, crossover_opts, update_opts):
    # expecting pmin, pmax, F and s
    pmin = (
        1 / pow(len(new_pop[0].x), 2)
        if "pmin" not in update_opts
        else update_opts["pmin"]
    )  # 1/n^2
    pmax = 1 / 2 if "pmax" not in update_opts else update_opts["pmax"]
    F = 1 + 0.2 if "F" not in update_opts else update_opts["F"]  # or 1 + 0.1
    s = math.e - 1 if "s" not in update_opts else update_opts["s"]  # e-1
    mutation_opts["p"] = (
        1 / len(new_pop[0].x) if "p" not in mutation_opts else mutation_opts["p"]
    )  # 1/n is default p0
    # https://arxiv.org/pdf/1902.02588
    # Only for 1+1 EA
    assert len(old_pop) == len(new_pop)
    assert len(old_pop) == 1
    if new_pop[0].fitness <= old_pop[0].fitness:
        mutation_opts["p"] = min(pow(F, s) * mutation_opts["p"], pmax)
    else:
        mutation_opts["p"] = max(mutation_opts["p"] / F, pmin)

    return mutation_opts, crossover_opts
