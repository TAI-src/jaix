from pymoo.util.ref_dirs import get_reference_directions
from joblib import Memory
from pymoo.core.problem import Problem as PyMOOProblem

memory = Memory("cache", verbose=0)


@memory.cache
def get_ref_dirs(m: int, num_refpoints: str | int):
    if isinstance(num_refpoints, int):
        n_refpoints = num_refpoints
    else:
        n_refpoints = get_num_refpoints(m, num_refpoints)
    ref_dirs = get_reference_directions("energy", m, n_refpoints)
    return ref_dirs


def get_num_refpoints(m: int, method: str):
    if method == "original":
        lookup_dict = {3: 91, 5: 210, 8: 156, 10: 275, 15: 135}
    elif method == "energy_small":
        lookup_dict = {3: 50, 5: 100, 8: 200, 10: 300, 15: 300}
    elif method == "energy_medium":
        lookup_dict = {3: 250, 5: 250, 8: 500, 10: 600, 15: 600}
    elif method == "energy_large":
        lookup_dict = {3: 500, 5: 500, 8: 1000, 10: 1000, 15: 1000}
    else:
        raise ValueError(f"Unsupported method: {method!r}")
    if m < min(lookup_dict.keys()) or m > max(lookup_dict.keys()):
        raise ValueError(f"Unsupported number of objectives: {m!r}")
    ref_points = lookup_dict.get(m, None)
    if ref_points is None:
        # interpolate between the two closest values
        keys = sorted(lookup_dict.keys())
        for i in range(len(keys) - 1):
            if keys[i] < m < keys[i + 1]:
                ref_points = int(
                    lookup_dict[keys[i]]
                    + (lookup_dict[keys[i + 1]] - lookup_dict[keys[i]])
                    * (m - keys[i])
                    / (keys[i + 1] - keys[i])
                )
                break
    if ref_points is None:
        raise ValueError(f"Unsupported number of objectives: {m!r}")
    return ref_points


@memory.cache
def get_pf(problem: PyMOOProblem, num_refpoints: str | int = "energy_large"):
    ref_dirs = get_ref_dirs(problem.n_obj, num_refpoints)
    try:
        pf = problem.pareto_front()
    except Exception as e:
        print(e)
        pf = problem.pareto_front(ref_dirs=ref_dirs)
    ideal = pf.min(axis=0)
    nadir = pf.max(axis=0) + 1e-6  # add small value to avoid zero range
    return pf, ideal, nadir
