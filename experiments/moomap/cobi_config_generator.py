## This file is adapted from the Cobi-problem-generator repository. The original file can be found here:
# https://github.com/numbbo/cobi-problem-generator/blob/6042a433387b36f718dc18108f69ba18aa863bb8/examples/simple_examples.py

import numpy as np
from jaix.env.utils.problem.cobi_problem import CobiProblemConfig
from ttex.config.config import ConfigFactory as CF


def _no_constraints():
    return {"Linear": [], "Quadratic": [], "Multi": []}


def _single_peak_objective(center, *, exponent=None, transformation=None):
    """Create one strictly convex-quadratic objective peak."""
    peak_transformation = None
    if exponent is not None:
        peak_transformation = {
            "name": "exponent",
            "params": {"exponent": exponent},
        }

    return {
        "H": np.array([np.eye(2)]),
        "c": np.array([center], dtype=float),
        "b": np.array([0.0]),
        "transformation": transformation,
        "peak_transformations": [peak_transformation],
    }


def _two_center_problem(
    *, peak_exponent=None, transformation_f1=None, transformation_f2=None
):
    """Problem used as the base geometry of Examples 1--4a."""
    objectives = (
        _single_peak_objective(
            [-2.0, 0.0],
            exponent=peak_exponent,
            transformation=transformation_f1,
        ),
        _single_peak_objective(
            [2.0, 0.0],
            exponent=peak_exponent,
            transformation=transformation_f2,
        ),
    )

    return {
        "n_var": 2,
        "objectives": objectives,
        "constraints": _no_constraints(),
        "domain": (-3.0, 3.0),
        "boundary_constraints": False,
    }


# -----------------------------------------------------------------------------
# 1. Linear Pareto front
# -----------------------------------------------------------------------------


def create_linear_front_problem():
    """The straight line Pareto front."""
    return _two_center_problem(peak_exponent=0.5)


# -----------------------------------------------------------------------------
# 2. Convex Pareto front
# -----------------------------------------------------------------------------


def create_convex_front_problem():
    """The convex front."""
    return _two_center_problem()


# -----------------------------------------------------------------------------
# 3. Concave Pareto front
# -----------------------------------------------------------------------------


def create_concave_front_problem():
    """The concave Pareto front."""
    return _two_center_problem(peak_exponent=0.25)


# -----------------------------------------------------------------------------
# 4a. Disconnected Pareto front using a transformation
# -----------------------------------------------------------------------------


def create_disconnected_linear_front_problem():
    """Two disconnected linear Pareto-front pieces created by a step transformation."""
    step = {
        "name": "step",
        "params": {
            "value": np.sqrt(2.0),
            "shift": 0.8,
        },
    }
    return _two_center_problem(
        peak_exponent=0.5,
        transformation_f1=step,
    )


# -----------------------------------------------------------------------------
# 4b. Disconnected Pareto front from multiple quadratic peaks
# -----------------------------------------------------------------------------


def create_disconnected_convex_parts_problem():
    """
    Two disconnected convex Pareto-front parts, using two quadratic peaks.
    """
    centers_f1 = np.array(
        [
            [-1.0, 3.0],
            [-1.0, -3.0],
        ]
    )
    centers_f2 = np.array(
        [
            [1.0, 3.0],
            [1.0, -3.0],
        ]
    )

    hessians = np.repeat(np.eye(2)[None, :, :], 2, axis=0)

    objectives = (
        {
            "H": hessians.copy(),
            "c": centers_f1,
            "b": np.array([0.0, 3.0]),
            "transformation": None,
            "peak_transformations": [None, None],
        },
        {
            "H": hessians.copy(),
            "c": centers_f2,
            "b": np.array([3.0, 0.0]),
            "transformation": None,
            "peak_transformations": [None, None],
        },
    )

    return {
        "n_var": 2,
        "objectives": objectives,
        "constraints": _no_constraints(),
        "domain": (-4.0, 4.0),
        "boundary_constraints": False,
    }


# -----------------------------------------------------------------------------
# 5. Many local Pareto fronts
# -----------------------------------------------------------------------------


def create_many_local_fronts_problem():
    """A simple 4-by-4 multipeak problem with 16 local peak-pair fronts."""
    y = np.array([-3.0, -1.0, 1.0, 3.0])

    centers_f1 = np.column_stack((-3.0 * np.ones(4), y))
    centers_f2 = np.column_stack((3.0 * np.ones(4), y))

    # Keeping every Hessian equal to I makes the construction easy to explain;
    # the different center separations and value shifts create distinct local
    # fronts in objective space.
    hessians = np.repeat(np.eye(2)[None, :, :], 4, axis=0)

    objectives = (
        {
            "H": hessians.copy(),
            "c": centers_f1,
            "b": np.array([0.0, 1.0, 2.0, 3.0]),
            "transformation": None,
            "peak_transformations": [None] * 4,
        },
        {
            "H": hessians.copy(),
            "c": centers_f2,
            "b": np.array([3.0, 2.0, 1.0, 0.0]),
            "transformation": None,
            "peak_transformations": [None] * 4,
        },
    )

    return {
        "n_var": 2,
        "objectives": objectives,
        "constraints": _no_constraints(),
        "domain": (-4.0, 4.0),
        "boundary_constraints": False,
    }


# -----------------------------------------------------------------------------
# 6. Few local Pareto fronts
# -----------------------------------------------------------------------------


def create_few_local_fronts_problem(peak_exponent=0.25):
    """A simple 2-by-2 multipeak problem with 4 local peak-pair fronts."""
    y = np.array([-3.0, 3.0])

    centers_f1 = np.column_stack((-3.0 * np.ones(2), y))
    centers_f2 = np.column_stack((3.0 * np.ones(2), y))

    # Keeping every Hessian equal to I makes the construction easy to explain;
    # the different center separations and value shifts create distinct local
    # fronts in objective space.
    hessians = np.repeat(np.eye(2)[None, :, :], 2, axis=0)

    peak_transformation = {
        "name": "exponent",
        "params": {"exponent": peak_exponent},
    }

    objectives = (
        {
            "H": hessians.copy(),
            "c": centers_f1,
            "b": np.array([0.0, 1.0]),
            "transformation": peak_transformation,
            "peak_transformations": [None] * 2,
        },
        {
            "H": hessians.copy(),
            "c": centers_f2,
            "b": np.array([1.0, 0.0]),
            "transformation": None,
            "peak_transformations": [None] * 2,
        },
    )

    return {
        "n_var": 2,
        "objectives": objectives,
        "constraints": _no_constraints(),
        "domain": (-4.0, 4.0),
        "boundary_constraints": False,
    }


def get_config(func_id: int):
    """Get the configuration for a specific example problem based on its function ID."""
    func_map = {
        1: create_linear_front_problem,
        2: create_convex_front_problem,
        3: create_concave_front_problem,
        4: create_disconnected_linear_front_problem,
        5: create_disconnected_convex_parts_problem,
        6: create_many_local_fronts_problem,
        7: create_few_local_fronts_problem,
    }
    if func_id not in func_map:
        raise ValueError(f"Invalid function ID {func_id}. Must be between 1 and 7.")
    config_dict = {
        "jaix.env.utils.problem.cobi_problem.CobiProblemConfig": func_map[func_id]()
    }

    return CF.from_dict(config_dict)


def get_configs():
    """Get configurations for all example problems."""
    return [get_config(func_id) for func_id in range(1, 8)]
