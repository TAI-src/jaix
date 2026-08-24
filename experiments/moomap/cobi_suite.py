from jaix.suite.ec_suite import ECSuite, ECSuiteConfig
from jaix.env.utils.problem.cobi_problem import CobiProblemConfig, CobiProblem


def get_cobi_config(type: str, n_var: int) -> CobiProblemConfig:
    if type == "moomap":
        return MoomapCobiProblemConfig(n_var=n_var)
    else:
        raise ValueError(f"Unknown cobi config type: {type}")


class MoomapCobiProblemConfig(CobiProblemConfig):
    def __init__(
        n_var: int = 10,
    ):
        super().__init__(
            n_var=n_var,
            domain=(-5.0, 5.0),
            n_peaks=10,
            peaks_value_shift=1.0,
            peaks_condition_number=100.0,
            peaks_alphas=[1.0] * 10,
            alpha=1.0,
            n_constraints={"Linear": 0, "Quadratic": 0, "Multi": 0},
            boundary_constraints=False,
            quadratic_constraints_size=1.0,
            quadratic_constraints_condition_number=None,
            n_multi_constraints_groups=0,
            n_multi_constraints_group_linear=0,
            n_multi_constraints_group_quadratic=0,
            constraints_feasible=True,
            perpendicular_linear_constraints=False,
            n_digits=None,
        )


# TODO: make a factory for creating the cobi problem configs
#
#
# TODO: Make a factory for creating a


# TODO: just turn this into a factory for creating the cobi problem configs, and then use that in the regular EC config
