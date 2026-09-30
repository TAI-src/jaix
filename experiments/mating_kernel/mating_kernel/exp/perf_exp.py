from pymoo.algorithms.moo.nsga2 import NSGA2
from ttex.config import Config, ConfigurableObject

from jaix.env.utils.problem.static_problem import StaticProblem


class PerfExperimentConfig(Config):
    def __init__(self, algorithm_name: str):
        super().__init__()
        self.algorithm_name = algorithm_name

    @staticmethod
    def init_algorithm(algorithm_name: str, problem: StaticProblem):

        if algorithm_name == "nsga2":

            return NSGA2

        else:
            raise ValueError(
                f"Invalid algorithm name: {algorithm_name}. Choose 'nsga2' or 'nsga3'."
            )


class PerfExperiment(ConfigurableObject):
    def __init__(self, config):
        super().__init__(config)
        self.config = config
