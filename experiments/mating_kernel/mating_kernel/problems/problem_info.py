from jaix.env.utils.problem.static_problem import StaticProblem


class ProblemInfo:
    def __init__(self, problem: StaticProblem):
        self.problem = str(problem)
        self.ideal_point = getattr(problem, "ideal_point", None)
        self.nadir_point = getattr(problem, "nadir_point", None)
        self.num_objectives = problem.num_objectives
        self.num_variables = problem.dimension
        self.lower_bounds = problem.lower_bounds
        self.upper_bounds = problem.upper_bounds
        if hasattr(
            problem, "name"
        ):  # CobiProblem has a name attribute, while REProblem does not
            self.uuid = problem.name
        elif hasattr(
            problem, "problem_name"
        ):  # REProblem has a problem_name attribute, while CobiProblem does not
            self.uuid = problem.problem_name
        else:  # Fallback to a generic UUID if neither attribute is present
            self.uuid = f"{self.problem}_{self.num_variables}_{self.num_objectives}"
            if hasattr(problem, "inst"):
                self.uuid += f"_{problem.inst}"
