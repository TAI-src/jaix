
from cobi import create_random_problem
from ttex.config import Config

from jaix.env.utils.problem.cobi_problem import CobiProblemConfig


class RandomCobiProblemConfig(Config):
    """
    Defaults are taken from https://github.com/numbbo/cobi-problem-generator/blob/6042a433387b36f718dc18108f69ba18aa863bb8/src/cobi/problem_generator.py#L208
    """

    def __init__(
        self,
        n_var: int = 2,
        domain: tuple[float, float] = (-5, 5),
        n_peaks: tuple[int | tuple[int, int], int | tuple[int, int]] = (
            (2, 5),
            (2, 5),
        ),
        peaks_value_shift: float | tuple[float, float] = 10,
        peaks_condition_number: float | tuple[float, float] | None = None,
        n_constraints: dict[str, int | tuple[int, int]] | None = None,
        boundary_constraints: bool = True,
        quadratic_constraints_size: float | tuple[float, float] = 10,
        quadratic_constraints_condition_number: float | tuple[float, float] | None = None,
        n_multi_constraints_groups: int | tuple[int, int] = 2,
        n_multi_constraints_group_linear: int | tuple[int, int] = (0, 1),
        n_multi_constraints_group_quadratic: int | tuple[int, int] = (2, 3),
        constraints_feasible: bool = True,
        perpendicular_linear_constraints: bool = False,
        n_digits: int | None = None,
    ):
        super().__init__()
        self.n_var = n_var
        self.domain = domain
        self.n_peaks = n_peaks
        self.peaks_value_shift = peaks_value_shift
        self.peaks_condition_number = peaks_condition_number

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


class CobiProblemConfigFactory:

    @staticmethod
    def create(config: RandomCobiProblemConfig, inst: int) -> CobiProblemConfig:
        random_cobi_problem = create_random_problem(
            n_var=config.n_var,
            seed=inst,
            domain=config.domain,
            n_peaks=config.n_peaks,
            peaks_value_shift=config.peaks_value_shift,
            peaks_condition_number=config.peaks_condition_number,
            n_constraints=config.n_constraints,
            boundary_constraints=config.boundary_constraints,
            quadratic_constraints_size=config.quadratic_constraints_size,
            quadratic_constraints_condition_number=config.quadratic_constraints_condition_number,
            n_multi_constraints_groups=config.n_multi_constraints_groups,
            n_multi_constraints_group_linear=config.n_multi_constraints_group_linear,
            n_multi_constraints_group_quadratic=config.n_multi_constraints_group_quadratic,
            constraints_feasible=config.constraints_feasible,
            perpendicular_linear_constraints=config.perpendicular_linear_constraints,
            n_digits=config.n_digits,
            print_seed=False,
        )
        config = CobiProblemConfig(
            n_var=random_cobi_problem.n_var,
            objectives=random_cobi_problem.objectives,
            constraints=random_cobi_problem.constraints,
            domain=random_cobi_problem.domain,
            alpha=(1, 1),  # Default value according to current cobi repo version
            boundary_constraints=random_cobi_problem.boundary_constraints,
        )
        return config
