def n_var(problem_name, n_obj, k_override: int | None = None):
    if k_override is not None:
        k = k_override
    elif problem_name == "dtlz1":
        k = 5
    elif problem_name == "dtlz7":
        k = 20
    elif problem_name in ["dtlz2", "dtlz3", "dtlz4", "dtlz5", "dtlz6"]:
        k = 10
    else:
        raise ValueError(f"Unsupported problem: {problem_name!r}")
    n_var = n_obj + k - 1
    return n_var
