import os
from typing import ClassVar

import numpy as np
from ttex.config import Config, ConfigurableObject

from jaix.env.utils.problem.re_problem.reproblem import (
    CRE21,
    CRE22,
    CRE23,
    CRE24,
    CRE25,
    CRE31,
    CRE32,
    CRE51,
    RE21,
    RE22,
    RE23,
    RE24,
    RE25,
    RE31,
    RE32,
    RE33,
    RE34,
    RE35,
    RE36,
    RE37,
    RE41,
    RE42,
    RE61,
    RE91,
)
from jaix.env.utils.problem.static_problem import StaticProblem


class REProblemConfig(Config):
    # We don't actually need a config for REProblem, but we need to define one to satisfy the ConfigurableObject interface
    def __init__(self):
        Config.__init__(self)


class REProblem(ConfigurableObject, StaticProblem):
    config_class = REProblemConfig
    problem_map: ClassVar[dict[str, type]] = {
        "RE21": RE21,  # 0
        "RE22": RE22,  # 1
        "RE23": RE23,  # 2
        "RE24": RE24,  # 3
        "RE25": RE25,  # 4
        "RE31": RE31,  # 5
        "RE32": RE32,  # 6
        "RE33": RE33,  # 7
        "RE34": RE34,  # 8
        "RE35": RE35,  # 9
        "RE36": RE36,  # 10
        "RE37": RE37,  # 11
        "RE41": RE41,  # 12
        "RE42": RE42,  # 13
        "RE61": RE61,  # 14
        "RE91": RE91,  # 15
        "CRE21": CRE21,  # 16
        "CRE22": CRE22,  # 17
        "CRE23": CRE23,  # 18
        "CRE24": CRE24,  # 19
        "CRE25": CRE25,  # 20
        "CRE31": CRE31,  # 21
        "CRE32": CRE32,  # 22
        "CRE51": CRE51,  # 23
    }
    unconstrained_version_map: ClassVar[dict[str, str]] = {
        "CRE21": "RE31",
        "CRE22": "RE32",
        "CRE23": "RE33",
        "CRE24": "RE34",
        "CRE25": "RE35",
        "CRE31": "RE41",
        "CRE32": "RE42",
        "CRE51": "RE61",
    }

    def __init__(self, config: REProblemConfig, inst: int):
        ConfigurableObject.__init__(self, config)
        if inst < 0 or inst >= len(REProblem.problem_map):
            raise ValueError(
                f"Invalid instance number {inst}. Must be between 0 and {len(REProblem.problem_map)-1}."
            )
        self.problem_name = list(REProblem.problem_map.keys())[inst]
        self.re_problem = REProblem.problem_map[
            self.problem_name
        ]()  # Instantiate the problem class
        self.constrained = (
            getattr(self.re_problem, "n_constraints", 0) > 0
        )  # Determine if the problem is constrained
        if self.constrained:
            self.unconstrained_problem_name = REProblem.unconstrained_version_map[
                self.problem_name
            ]
        else:
            self.unconstrained_problem_name = self.problem_name

        # Initialize StaticProblem accordingly
        self.lower_bounds = getattr(self.re_problem, "lbound", [])
        assert len(self.lower_bounds) > 0
        self.upper_bounds = getattr(self.re_problem, "ubound", [])
        StaticProblem.__init__(
            self,
            dimension=getattr(self.re_problem, "n_variables", -1),
            num_objectives=getattr(self.re_problem, "n_objectives", -1),
            precision=None,
        )

        # Retrieve the ideal and nadir point from the files.
        self._ideal_point, self._nadir_point = self.extract_ideal_nadir()

    @property
    def nadir_point(self) -> np.ndarray:
        return self._nadir_point

    @property
    def ideal_point(self) -> np.ndarray:
        return self._ideal_point

    def extract_ideal_nadir(self):
        # get the absolute path of this file so the data files can be relative
        this_folder = os.path.dirname(os.path.abspath(__file__))
        folder = os.path.join(this_folder, "ideal_nadir_points")

        ideal_file = os.path.join(
            folder, f"ideal_point_{self.unconstrained_problem_name}.dat"
        )
        nadir_file = os.path.join(
            folder, f"nadir_point_{self.unconstrained_problem_name}.dat"
        )

        ideal_point = np.loadtxt(ideal_file)
        nadir_point = np.loadtxt(nadir_file)

        # For the unconstrained version, the last objective is the constraint violation, which we don't want to include in the ideal and nadir points. So we remove it.
        if self.constrained:
            ideal_point = ideal_point[:-1]
            nadir_point = nadir_point[:-1]

        return ideal_point, nadir_point

    def _eval(self, x):
        if self.constrained:
            f, g = self.re_problem.evaluate(x)
            if sum(g) >= 0:  # There is a constraint violation
                f = [np.nan] * self.num_objectives
        else:
            f = self.re_problem.evaluate(x)
        return f, f  # Return the same values since all of these are not noisy

    def __str__(self) -> str:
        return f"REProblem_{self.problem_name}"
