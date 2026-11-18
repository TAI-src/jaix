from typing import Self

from pydantic import BaseModel

from jaix.env.utils.problem.static_problem import StaticProblem


class ProblemInfo(BaseModel):
    problem: str
    ideal_point: list[float] | None = None
    nadir_point: list[float] | None = None
    num_objectives: int
    num_variables: int
    lower_bounds: list[float]
    upper_bounds: list[float]
    uuid: str

    @classmethod
    def from_problem(cls, problem: StaticProblem) -> Self:
        num_variables = problem.dimension
        num_objectives = problem.num_objectives

        uuid = getattr(problem, "name", None) or getattr(problem, "problem_name", None)

        if uuid is None:
            uuid = f"{problem}_{num_variables}_{num_objectives}"
            if hasattr(problem, "inst"):
                uuid += f"_{problem.inst}"

        return cls(
            problem=str(problem),
            ideal_point=(
                list(getattr(problem, "ideal_point"))  # noqa: B009
                if getattr(problem, "ideal_point", None) is not None
                else None
            ),
            nadir_point=(
                list(getattr(problem, "nadir_point"))  # noqa: B009
                if getattr(problem, "nadir_point", None) is not None
                else None
            ),
            num_objectives=num_objectives,
            num_variables=num_variables,
            lower_bounds=list(problem.lower_bounds),
            upper_bounds=list(problem.upper_bounds),
            uuid=uuid,
        )
