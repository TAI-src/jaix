import numpy as np
from jaix.env.utils.problem.static_problem import StaticProblem
from pymoo.core.problem import ElementwiseProblem


class PymooProblemWrapper(ElementwiseProblem):
    def __init__(self, static_problem: StaticProblem, record: bool = True):
        self.static_problem = static_problem
        super().__init__(
            n_var=static_problem.dimension,
            n_obj=static_problem.num_objectives,
            n_ieq_constr=0,
            xl=static_problem.lower_bounds,
            xu=static_problem.upper_bounds,
        )
        self.records: list[dict] = []
        self.record = record

    def _evaluate(self, X, out, *args, **kwargs):
        F, _ = self.static_problem(X)
        out["F"] = np.array(F)
        if self.record:
            self.records.append({"X": X, "F": F})

    def retrieve_records(self):
        records = self.records
        self.records = []
        return records
