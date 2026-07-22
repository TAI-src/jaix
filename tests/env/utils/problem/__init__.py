try:
    from jaix.env.utils.problem.cobi_problem import CobiProblem, CobiProblemConfig
except ImportError:
    # If the import fails, we set CobiProblem and CobiProblemConfig to None
    CobiProblem = None  # type: ignore[assignment,misc]
    CobiProblemConfig = None  # type: ignore[assignment,misc]
