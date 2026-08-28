try:
    from jaix.env.utils.problem.cobi_problem import CobiProblem, CobiProblemConfig
    from jaix.env.utils.problem.cobi_problem_factory import (
        CobiProblemConfigFactory,
        RandomCobiProblemConfig,
    )
except ImportError:
    # If the import fails, we set CobiProblem and CobiProblemConfig to None
    CobiProblem = None  # type: ignore[assignment,misc]
    CobiProblemConfig = None  # type: ignore[assignment,misc]
    CobiProblemConfigFactory = None  # type: ignore[assignment,misc]
    RandomCobiProblemConfig = None  # type: ignore[assignment,misc]
