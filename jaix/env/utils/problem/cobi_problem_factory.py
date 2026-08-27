from typing import Dict, Optional, Tuple, Union
from ttex.config import Config
from cobi import create_random_problem
from jaix.env.utils.problem.cobi_problem import CobiProblemConfig, CobiProblem


class RandomCobiProblemConfig(Config):
    """
    Defaults are taken from https://github.com/numbbo/cobi-problem-generator/blob/6042a433387b36f718dc18108f69ba18aa863bb8/src/cobi/problem_generator.py#L208
    """

    def __init__(
        self,
        n_var: int = 2,
        domain: Tuple[float, float] = (-5, 5),
        n_peaks: Tuple[Union[int, Tuple[int, int]], Union[int, Tuple[int, int]]] = (
            (2, 5),
            (2, 5),
        ),
        peaks_value_shift: Union[float, Tuple[float, float]] = 10,
        peaks_condition_number: Optional[Union[float, Tuple[float, float]]] = None,
        n_constraints: Optional[Dict[str, Union[int, Tuple[int, int]]]] = None,
        boundary_constraints: bool = True,
        quadratic_constraints_size: Union[float, Tuple[float, float]] = 10,
        quadratic_constraints_condition_number: Optional[
            Union[float, Tuple[float, float]]
        ] = None,
        n_multi_constraints_groups: Union[int, Tuple[int, int]] = 2,
        n_multi_constraints_group_linear: Union[int, Tuple[int, int]] = (0, 1),
        n_multi_constraints_group_quadratic: Union[int, Tuple[int, int]] = (2, 3),
        constraints_feasible: bool = True,
        perpendicular_linear_constraints: bool = False,
        n_digits: Optional[int] = None,
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
