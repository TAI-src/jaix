from pymoo.core.population import Population
import numpy as np
from pymoo.problems import get_problem


def dummy_comp(pop, P, **kwargs):
    return P[:, 0]


def create_pop(problem=None, size=10):
    if problem is None:
        problem = get_problem("zdt1")
    if hasattr(problem, "xl") and hasattr(problem, "xu"):
        x = np.random.uniform(
            low=problem.xl, high=problem.xu, size=(size, len(problem.xl))
        )
    elif hasattr(problem, "lower_bounds") and hasattr(problem, "upper_bounds"):
        x = np.random.uniform(
            low=problem.lower_bounds,
            high=problem.upper_bounds,
            size=(size, len(problem.lower_bounds)),
        )
    else:
        raise ValueError("Problem does not have defined lower and upper bounds.")
    pop = Population.new("X", x)
    if hasattr(problem, "evaluate"):
        pop.set("F", problem.evaluate(pop.get("X")))
    else:
        for ind in pop:
            _, y = problem(ind.X)
            ind.F = y
    return pop
