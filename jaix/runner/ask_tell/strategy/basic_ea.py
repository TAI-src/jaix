from ttex.config import Config, ConfigurableObject
from enum import Enum
from jaix.runner.ask_tell import ATStrategy

EAStrategy = Enum("EAStrategy", [("Comma", 0), ("Plus", 1)])


class MutationOp(Enum):
    Gaussian = 0
    Uniform = 1
    Polynomial = 2
    BitFlip = 3
    Permutation = 4
    Swap = 5
    Inversion = 6
    Scramble = 7
    Heuristic = 8
    Random = 9


class CrossoverOp(Enum):
    OnePoint = 0
    TwoPoint = 1
    Uniform = 2
    Heuristic = 3
    Random = 4


class BasicEAConfig(Config):
    def __init__(
        self,
        strategy: EAStrategy,
        population_size: int,
        num_parents: int,
        mutation_op: MutationOp,
        crossover_op: CrossoverOp,
        mutation_opts={},
        crossover_opts={},
    ):
        self.strategy = strategy
        self.population_size = population_size
        self.num_parents = num_parents
        self.mutation_op = mutation_op
        self.mutation_opts = mutation_opts
        self.crossover_op = crossover_op
        self.crossover_opts = crossover_opts


class BasicEA(ConfigurableObject, ATStrategy):
    def __init__(self, config: BasicEAConfig, xstart, *args, **kwargs):
        ConfigurableObject.__init__(self, config)
        ATStrategy.__init__(self, xstart)
