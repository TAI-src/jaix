from jaix.env.utils.problem.static_problem import StaticProblem


class ProblemInfo:
    def __init__(self, problem: StaticProblem):
        self.problem = str(problem)
        self.ideal_point = problem.ideal_point
        self.nadir_point = problem.nadir_point
        self.num_objectives = problem.num_objectives
        self.num_variables = problem.dimension
        self.lower_bounds = problem.lower_bounds
        self.upper_bounds = problem.upper_bounds
        self.problem_id = f"{self.problem}_{problem.dimension}_{problem.num_objectives}_{problem.inst}"
