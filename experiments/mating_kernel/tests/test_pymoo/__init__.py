from pymoo.core.population import Population
import numpy as np
from pymoo.problems import get_problem


def dummy_comp(pop, P, **kwargs):
    return P[:, 0]


def create_pop(problem=None, size=10):
    if problem is None:
        problem = get_problem("zdt1")
    x = np.random.uniform(low=problem.xl, high=problem.xu, size=(size, len(problem.xl)))
    pop = Population.new("X", x)
    pop.set("F", problem.evaluate(pop.get("X")))
    return pop
