from jaix.env.utils.problem.static_problem import StaticProblem
from pymoo.core.problem import ElementwiseProblem
import numpy as np
from jaix.env.utils.archive.mo_archive import MOArchiveConfig, MOArchive, KeepDominated
from config_nsga3x import MOEvalEntry
from jaix.env.utils.archive.entry_scorer import (
    ReferenceVectorDistanceScorer,
)
from jaix.env.singular.ec_env import ECEnvironment, ECEnvironmentConfig


class PymooProblemWrapper(ElementwiseProblem):
    def __init__(self, static_problem: StaticProblem):
        self.static_problem = static_problem
        super().__init__(
            n_var=static_problem.dimension,
            n_obj=static_problem.num_objectives,
            n_ieq_constr=0,
            xl=static_problem.lower_bounds,
            xu=static_problem.upper_bounds,
        )
        self.archive = PymooProblemWrapper._create_eval_archive(static_problem)

    @staticmethod
    def _create_eval_archive(func: StaticProblem) -> MOArchive:
        archive_config = MOArchiveConfig(
            MOEvalEntry,
            secondary_criterion_class=ReferenceVectorDistanceScorer,
            max_size=None,
            keep_dominated=KeepDominated.NONE,
            only_new_entries=False,
            num_refpoints="original",
        )
        env = ECEnvironment(ECEnvironmentConfig(budget_multiplier=1), func=func)
        return MOArchive(archive_config, env=env)

    def _evaluate(self, X, out, *args, **kwargs):
        F, _ = self.static_problem(X)
        entry = MOEvalEntry(x=np.array(X), y=np.array(F))
        self.archive.add([entry])
        out["F"] = np.array(F)
