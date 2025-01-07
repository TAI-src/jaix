from ttex.config import Config, ConfigurableObject
from enum import Enum
from jaix.runner.ask_tell import ATStrategy
import numpy as np
from jaix.runner.ask_tell.strategy.utils.ea_utils import (
    global_flip,
    onepoint_crossover,
    uniform_crossover,
    Individual,
    select,
)
from jaix.env.composite import CompositeEnvironment
from typing import Optional

EAStrategy = Enum("EAStrategy", [("Comma", 0), ("Plus", 1)])


class MutationOp(Enum):
    FLIP = global_flip


class CrossoverOp(Enum):
    ONEPOINT = onepoint_crossover
    UNIFORM = uniform_crossover


class BasicEAConfig(Config):
    def __init__(
        self,
        strategy: EAStrategy,
        mu: int,
        lam: int,
        mutation_op: Optional[MutationOp],
        crossover_op: Optional[CrossoverOp],
        mutation_opts={},
        crossover_opts={},
    ):
        self.strategy = strategy
        self.mu = mu
        self.lam = lam
        self.mutation_op = mutation_op
        self.mutation_opts = mutation_opts
        self.crossover_op = crossover_op
        self.crossover_opts = crossover_opts
        if self.strategy == EAStrategy.Comma:
            assert self.lam >= self.mu
        if self.crossover_op is None:
            assert self.lam <= self.mu


class BasicEA(ConfigurableObject, ATStrategy):
    def __init__(self, config: BasicEAConfig, xstart, *args, **kwargs):
        ConfigurableObject.__init__(self, config)
        print(len(xstart))
        print(xstart)
        assert len(xstart) == config.mu
        ATStrategy.__init__(self, xstart)

    def initialize(self):
        self.pop = [Individual(np.array(xi), np.nan) for xi in self.xstart]

    @property
    def name(self):
        return f"Basic EA (mu={self.mu} {self.strategy} lam={self.lam})"

    def ask(self, env, **kwargs):
        offspring = [None] * self.lam
        # Variation
        for i in range(self.lam):
            if self.crossover_op is not None:
                parents_idx = np.random.randint(low=0, high=self.mu, size=2)
                child_x = self.crossover_op(
                    self.pop[parents_idx[0]].x,
                    self.pop[parents_idx[1]].x,
                    **self.crossover_opts,
                )
            else:
                child_x = self.pop[i].x
            if self.mutation_op is not None:
                child_x = self.mutation_op(child_x, **self.mutation_opts)
            offspring[i] = child_x
        return offspring

    def tell(self, env, solutions, function_values, **kwargs):
        assert len(solutions) == len(function_values)
        if isinstance(env.unwrapped, CompositeEnvironment):
            function_values = [v for n, v in function_values]
        # TODO: currently only doing single-objective
        # TODO: make this setup common to avoid code duplication
        assert all([len(v) == 1 for v in function_values])
        new_pop = [Individual(x, f[0]) for x, f in zip(solutions, function_values)]
        # Survival selection
        if self.strategy == EAStrategy.Comma:
            self.pop = select(new_pop, self.mu)
        elif self.strategy == EAStrategy.Plus:
            self.pop += new_pop
            self.pop = select(self.pop + new_pop, self.mu)
