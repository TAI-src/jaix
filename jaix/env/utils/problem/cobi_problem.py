import numpy as np
from typing import Tuple, Union
from cobi import CobiProblem as OrigCobiProblem
from ttex.config import Config, ConfigurableObject

from jaix.env.utils.problem.static_problem import StaticProblem


class CobiProblemConfig(Config):
    """
    Defaults are taken from https://github.com/numbbo/cobi-problem-generator/blob/6042a433387b36f718dc18108f69ba18aa863bb8/src/cobi/cobi_problem.py#L449
    """

    def __init__(
        self,
        n_var: int,
        objectives: Tuple[dict, dict],
        constraints: dict[str, list] | None = None,
        domain: Tuple[float, float] = (-5, 5),
        alpha: Union[float, Tuple[float, float]] = (1, 1),
        boundary_constraints: bool = True,
    ):
        super().__init__()
        self.n_var = n_var
        self.objectives = objectives
        self.constraints = constraints
        self.domain = domain
        self.alpha = alpha
        self.boundary_constraints = boundary_constraints

        # convert info to StaticProblem attributes where possible
        self.lower_bounds = [self.domain[0]] * self.n_var
        self.upper_bounds = [self.domain[1]] * self.n_var


class CobiProblem(ConfigurableObject, StaticProblem):
    config_class = CobiProblemConfig

    def __init__(self, config: CobiProblemConfig, inst: int):
        ConfigurableObject.__init__(self, config)
        self.inst = inst
        self.cobi_problem = OrigCobiProblem(
            n_var=self.n_var,
            objectives=self.objectives,
            constraints=self.constraints,
            domain=self.domain,
            alpha=self.alpha,
            boundary_constraints=self.boundary_constraints,
        )

        StaticProblem.__init__(
            self,
            dimension=self.n_var,
            num_objectives=2,
        )

    @property
    def nadir_point(self) -> np.ndarray:
        return self.cobi_problem.nadir_point()

    @property
    def ideal_point(self) -> np.ndarray:
        return self.cobi_problem.ideal_point()

    def _eval(self, x):
        fitness = self.cobi_problem.evaluate_objectives(x)
        constraint_violations = self.cobi_problem.violation_point(x)

        if constraint_violations > 0:
            # Some constraint is violated
            return [np.nan] * self.num_objectives, [np.nan] * self.num_objectives
        else:
            return list(fitness), list(
                fitness
            )  # return both clean and noisy fitness (same values, cobi-problem is deterministic)
