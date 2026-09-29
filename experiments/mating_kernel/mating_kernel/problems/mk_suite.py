from ttex.config import Config, ConfigurableObject

from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from jaix.env.utils.problem.static_problem import StaticProblem
from mating_kernel.problems.cobi_configs import get_config
from mating_kernel.problems.cobi_configs import names as cobi_names
from mating_kernel.problems.problem_info import ProblemInfo


class MKSuiteConfig(Config):
    def __init__(
        self,
        cobi: bool = True,
        re: bool = True,
        num_objectives: list[int] | None = None,
        constrained: bool = False,
    ):
        super().__init__()
        self.cobi = cobi
        self.re = re
        self.num_objectives = num_objectives
        self.constrained = constrained


class MKSuite(ConfigurableObject):
    config_class = MKSuiteConfig

    def __init__(self, config: MKSuiteConfig):
        super().__init__(config)
        self.problems = MKSuite.generate_problem_list(
            self.cobi, self.re, self.constrained, self.num_objectives
        )
        self.problem_info_list = [ProblemInfo(p) for p in self.problems]
        self.problem_id_map = {
            pinfo.uuid: {"problem": p, "info": pinfo}
            for p, pinfo in zip(self.problems, self.problem_info_list)
        }

    @staticmethod
    def re_problem_list(constrained: bool = False) -> list[REProblem]:
        problem_max_id = 24 if constrained else 16
        # All non-constrained RE problems are included in the list.
        problems = [REProblem(REProblemConfig(), i) for i in range(problem_max_id)]
        return problems

    @staticmethod
    def cobi_problem_list() -> list[CobiProblem]:
        problem_max_id = 7
        cobi_configs = [get_config(func_id) for func_id in range(problem_max_id)]
        problems = [
            CobiProblem(config, inst=i) for i, config in enumerate(cobi_configs)
        ]
        for p in problems:
            p.name = cobi_names[p.inst]
        return problems

    @staticmethod
    def generate_problem_list(
        cobi: bool = True,
        re: bool = True,
        constrained: bool = False,
        num_objectives: list[int] | None = None,
    ) -> list[StaticProblem]:
        problems: list[StaticProblem] = []
        if cobi:
            cobi_problems = MKSuite.cobi_problem_list()
            problems.extend(cobi_problems)
        if re:
            re_problems = MKSuite.re_problem_list(constrained=constrained)
            problems.extend(re_problems)
        if num_objectives is not None:
            problems = [p for p in problems if p.num_objectives in num_objectives]
        return problems
