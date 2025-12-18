"""
A COCO problem from numbbo coco
"""

from cocoex import Problem
from jaix.env.utils.problem.static_problem import StaticProblem


class COCOProblem(StaticProblem):
    def __init__(self, problem: Problem):
        self.problem = problem
        StaticProblem.__init__(self, problem.dimension, problem.number_of_objectives)

    def __getattr__(self, name):
        """
        defer to the problems's attribute if not found here
        """
        if name != "problem":
            return self.problem.__getattribute__(name)
        else:
            return StaticProblem._getattr__(name)

    def __hasattr__(self, name):
        """
        defer to the problems's attributes if not found here
        """
        if name != "problem":
            return self.problem.__hasattribute__(name)
        else:
            return StaticProblem._hasattr__(name)

    def evalsleft(self, budget_multiplier):
        evals_allowed = self.dimension * budget_multiplier
        evals_done = max(self.evaluations, self.evaluations_constraints)
        evals_left = evals_allowed - evals_done
        return int(evals_left)

    def final_target_hit(self):
        return self.problem.final_target_hit

    def _eval(self, x):
        res = self.problem(x)
        if self.problem.number_of_objectives == 1:
            return [res], [res]
        else:
            return res, res

    def __str__(self):
        return f"{self.dimension}/{self.problem.name}"

    def close(self):
        self.problem.free()
        return super().close()
