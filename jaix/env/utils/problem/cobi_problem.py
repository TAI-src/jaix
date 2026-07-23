import numpy as np
from cobi import create_random_problem
from ttex.config import Config, ConfigurableObject

from jaix.env.utils.problem.static_problem import StaticProblem


class CobiProblemConfig(Config):
    """
    Defaults are taken from https://github.com/numbbo/cobi-problem-generator/blob/main/src/cobi/problem_generator.py
    """

    def __init__(
        self,
        n_var: int,  #  Dimension of the decision/search space.
        domain: tuple[float, float] = (
            -5,
            5,
        ),  # Lower and upper bounds for all decision variables.
        n_peaks: tuple[int | tuple[int, int], int | tuple[int, int]] = (
            (2, 5),
            (2, 5),
        ),  # Number of peaks for each objective function. If a tuple (min, max),
        # a random integer in that range is chosen.
        peaks_value_shift: float
        | tuple[
            float, float
        ] = 10,  #  Range from which the f-value shifts b1_i, b2_j for peaks are sampled.
        # If a single number x, shifts are sampled uniformly from [-x, x]. If a tuple (min, max), shifts are sampled uniformly from [min, max].
        peaks_condition_number: float
        | tuple[float, float]
        | None = None,  #  Condition number for Hessian matrices of objective functions.
        # If a number x, the actual condition number is sampled logarithmically from [1, x]. If a tuple (min, max), it is sampled from [min, max]. If None, Hessians with random condition numbers are generated.
        peaks_alphas: float
        | tuple[
            float, float
        ] = 1,  # Range from which the alphas for peaks are sampled. If a single number x,
        # all alphas are x. If a tuple (min, max), alphas are sampled uniformly from [min, max].
        alpha: float
        | tuple[float, float] = (
            1,
            1,
        ),  # Exponents used to transform the objective functions.
        # If a single number x, alpha_1 = alpha_2 = x.
        n_constraints: dict[str, int | tuple[int, int]]
        | None = None,  # Dictionary with the number of constraints of each type.
        # Must have the form {'Linear': number of linear constraints,
        # 'Quadratic': number of quadratic constraints,
        # 'Multi': number of multi-constraints}.
        boundary_constraints: bool = True,  # if True, automatically adds boundary constraints
        # for each decision variable, ensuring that constraint violations reflect the domain.
        quadratic_constraints_size: float
        | tuple[
            float, float
        ] = 10,  #  Range from which sizes of quadratic constraints are sampled.
        # If a single number x, all sizes are x. If a tuple (min, max),
        # sizes are sampled logarithmically from [min, max].
        quadratic_constraints_condition_number: float
        | tuple[float, float]
        | None = None,  # Condition number for Hessian matrices of quadratic constraints.
        # Follows the same rules as peaks_condition_number.
        n_multi_constraints_groups: int
        | tuple[int, int] = 2,  #  Number of groups in each multi-constraint.
        # A multi-constraint has the form min_k [max_l [g_{k,l}]] <= 0, where g_{k,l} are
        # linear or convex-quadratic constraints and k runs from 1 to the number of groups.
        # If a tuple (min, max) is provided, a random integer in that range is chosen
        # for each multi-constraint.
        n_multi_constraints_group_linear: int
        | tuple[int, int] = (
            0,
            1,
        ),  # Number of linear constraints in each group of a multi-constraint.
        # For group k of a multi-constraint, this defines the number
        # of linear constraints in {g_{k,1}, ..., g_{k,m_k}}.
        # If a tuple (min, max) is provided, a random integer in that range is chosen for each group.
        n_multi_constraints_group_quadratic: int
        | tuple[int, int] = (
            2,
            3,
        ),  # Number of quadratic constraints in each group of a multi-constraint.
        # For group k of a multi-constraint, this defines the number of
        # quadratic constraints in {g_{k,1}, ..., g_{k,m_k}}.
        # If a tuple (min, max) is provided, a random integer in that range is chosen for each group.
        constraints_feasible: bool = True,  # If True, all constraints are generated so that
        # some randomly sampled point is feasible.
        perpendicular_linear_constraints: bool = False,  # If True, linear constraints are
        # generated perpendicular to the x1-x2 plane.
        n_digits: int
        | None = None,  # If not None, rounds all generated numbers to this number of digits.
    ):
        Config.__init__(self)
        self.n_var = n_var
        self.domain = domain
        self.n_peaks = n_peaks
        self.peaks_value_shift = peaks_value_shift
        self.peaks_condition_number = peaks_condition_number
        self.peaks_alphas = peaks_alphas
        self.alpha = alpha
        self.n_constraints = n_constraints
        self.boundary_constraints = boundary_constraints
        self.quadratic_constraints_size = quadratic_constraints_size
        self.quadratic_constraints_condition_number = (
            quadratic_constraints_condition_number
        )
        self.n_multi_constraints_groups = n_multi_constraints_groups
        self.n_multi_constraints_group_linear = n_multi_constraints_group_linear
        self.n_multi_constraints_group_quadratic = n_multi_constraints_group_quadratic
        self.constraints_feasible = constraints_feasible
        self.perpendicular_linear_constraints = perpendicular_linear_constraints
        self.n_digits = n_digits

        # convert info to StaticProblem attributes where possible
        self.lower_bounds = [self.domain[0]] * self.n_var
        self.upper_bounds = [self.domain[1]] * self.n_var


class CobiProblem(ConfigurableObject, StaticProblem):
    config_class = CobiProblemConfig

    def __init__(self, config: CobiProblemConfig, inst: int):
        ConfigurableObject.__init__(self, config)
        self.inst = inst
        self.cobi_problem = create_random_problem(
            n_var=self.n_var,
            seed=self.inst,
            domain=self.domain,
            n_peaks=self.n_peaks,
            peaks_value_shift=self.peaks_value_shift,
            peaks_condition_number=self.peaks_condition_number,
            peaks_alphas=self.peaks_alphas,
            alpha=self.alpha,
            n_constraints=self.n_constraints,
            boundary_constraints=self.boundary_constraints,
            quadratic_constraints_size=self.quadratic_constraints_size,
            quadratic_constraints_condition_number=self.quadratic_constraints_condition_number,
            n_multi_constraints_groups=self.n_multi_constraints_groups,
            n_multi_constraints_group_linear=self.n_multi_constraints_group_linear,
            n_multi_constraints_group_quadratic=self.n_multi_constraints_group_quadratic,
            constraints_feasible=self.constraints_feasible,
            perpendicular_linear_constraints=self.perpendicular_linear_constraints,
            n_digits=self.n_digits,
            print_seed=False,  # don't print seed to avoid cluttering output
        )
        if self.n_digits is not None:
            # n_digits is the value the problem gets rounded to. So we should adapt precision accordingly
            precision = 10 ** (-self.n_digits)
        else:
            precision = None
        StaticProblem.__init__(
            self,
            dimension=self.n_var,
            num_objectives=2,
            precision=precision,
        )

    def _eval(self, x):
        fitness = self.cobi_problem.evaluate_objectives(x)
        constraint_violations = self.cobi_problem.violation_point(x)

        if constraint_violations > 0:
            # Some constraint is violated
            return [np.nan] * self.num_objectives, [np.nan] * self.num_objectives
        else:
            return list(fitness), list(
                fitness
            )  # return both clean and noisy fitness (same values, cobi-problem is deterministic)
