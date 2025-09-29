from ttex.config import Config, ConfigurableObject
from enum import Enum
from jaix.runner.ask_tell.at_strategy import ATStrategy
import numpy as np
from jaix.runner.ask_tell.strategy.utils.ea_utils import (
    global_flip,
    onepoint_crossover,
    uniform_crossover,
    Individual,
    select,
    ddl_update,
)
from jaix.env.composite.composite_environment import CompositeEnvironment
from typing import Optional
from gymnasium import Env, spaces

EAStrategy = Enum("EAStrategy", [("Comma", 0), ("Plus", 1)])

# Doing function strings instead of directly to ensure being able to configure via dict


class MutationOp(Enum):
    """
    Enum for mutation operators
    """

    FLIP = "global_flip"


class CrossoverOp(Enum):
    """
    Enum for crossover operators
    """

    ONEPOINT = "onepoint_crossover"
    UNIFORM = "uniform_crossover"


class UpdateStrategy(Enum):
    """
    Enum for update strategies
    """

    DDL = "ddl_update"


class WarmStartStrategy(Enum):
    """
    Enum for warm start strategies
    """

    BEST = "best"
    LAST = "last"
    NONE = None


class BasicEAConfig(Config):
    """
    Configuration class for BasicEA
    """

    def __init__(
        self,
        strategy: EAStrategy,  # Comma or Plus
        mu: int,  # number of parents
        lam: int,  # number of offspring
        mutation_op: Optional[MutationOp],  # mutation operator
        crossover_op: Optional[CrossoverOp],  # crossover operator
        mutation_opts={},  # mutation operator options
        crossover_opts={},  # crossover operator options
        warm_start_strategy: WarmStartStrategy = WarmStartStrategy.NONE,  # warm start
        update_strategy: Optional[UpdateStrategy] = None,  # update strategy
        update_opts={},  # update strategy options
    ):
        self.strategy = strategy
        self.mu = mu
        self.lam = lam
        self.mutation_op = mutation_op
        self.mutation_opts = mutation_opts
        self.crossover_op = crossover_op
        self.crossover_opts = crossover_opts
        self.warm_start_strategy = warm_start_strategy
        if self.strategy == EAStrategy.Comma:
            assert self.lam >= self.mu
        if self.crossover_op is None:
            assert self.lam <= self.mu
        else:
            assert self.mu >= 2
        assert self.mutation_op is not None or self.crossover_op is not None
        self.update_strategy = update_strategy
        self.update_opts = update_opts


class BasicEA(ConfigurableObject, ATStrategy):
    config_class = BasicEAConfig
    """
    Basic Evolutionary Algorithm
    """

    def __init__(self, config: BasicEAConfig, env: Env, *args, **kwargs):
        ConfigurableObject.__init__(self, config)
        self.xstart = [env.unwrapped.action_space.sample() for _ in range(self.mu)]
        # TODO: for now only MultiBinary or MultiDiscrete
        self.mutation_opts["low"] = [0] * len(self.xstart[0])
        if isinstance(env.unwrapped.action_space, spaces.MultiBinary):
            self.mutation_opts["high"] = [1] * len(self.xstart[0])
        elif isinstance(env.unwrapped.action_space, spaces.MultiDiscrete):
            # An nvec of 1 means only 1 option, so an empty dimension
            assert all(env.unwrapped.action_space.nvec > 1)
            self.mutation_opts["high"] = env.unwrapped.action_space.nvec - 1
        else:
            raise NotImplementedError(
                f"Action space {env.unwrapped.action_space} not supported yet"
            )
        self.mutate = (
            globals().get(self.mutation_op.value)
            if self.mutation_op is not None
            else None
        )
        self.crossover = (
            globals().get(self.crossover_op.value)
            if self.crossover_op is not None
            else None
        )
        self.update = (
            globals().get(self.update_strategy.value)
            if self.update_strategy is not None
            else None
        )

        ATStrategy.__init__(self, self.xstart)
        self.initialize()

    def initialize(self):
        """
        (Re-)Initialize the EA
        """
        self.gen = 0
        self.pop = [None] * self.mu
        assert len(self.xstart) == self.mu
        for i in range(self.mu):
            if not isinstance(self.xstart[i], Individual):
                self.pop[i] = Individual(np.array(self.xstart[i]), np.inf, self.gen)
            else:
                self.pop[i] = self.xstart[i]

    @property
    def name(self):
        return f"Basic EA (mu={self.mu} {self.strategy} lam={self.lam})"

    def ask(self, env, **kwargs):
        """
        Generate lam offspring
        """
        offspring = [None] * self.lam
        # Variation
        for i in range(self.lam):
            if self.crossover_op is not None:
                parents_idx = np.random.choice(
                    list(range(self.mu)), size=2, replace=False
                )
                child_x = self.crossover(
                    self.pop[parents_idx[0]].x,
                    self.pop[parents_idx[1]].x,
                    **self.crossover_opts,
                )
            else:
                child_x = self.pop[i].x
            if self.mutation_op is not None:
                child_x = self.mutate(child_x, **self.mutation_opts)
            offspring[i] = child_x
        return offspring

    def tell(self, env, solutions, function_values, **kwargs):
        """
        Update the EA with the new solutions and function values
        """
        assert len(solutions) == len(function_values)
        if isinstance(env.unwrapped, CompositeEnvironment):
            function_values = [v for n, v in function_values]
        # TODO: currently only doing single-objective
        # TODO: make this setup common to avoid code duplication
        assert all([len(v) == 1 for v in function_values])
        self.gen += 1
        new_pop = [
            Individual(x, f[0], self.gen) for x, f in zip(solutions, function_values)
        ]

        # Update parameters
        if self.update is not None:
            self.mutation_opts, self.crossover_opts = self.update(
                self.pop,
                new_pop,
                self.mutation_opts,
                self.crossover_opts,
                self.update_opts,
            )

        # Survival selection
        if self.strategy == EAStrategy.Comma:
            self.pop = select(new_pop, self.mu)
        elif self.strategy == EAStrategy.Plus:
            self.pop += new_pop
            self.pop = select(self.pop + new_pop, self.mu)

    def warm_start(self, xlast, env, **kwargs):
        """
        Warm start the strategy
        """
        self.xstart = [env.unwrapped.action_space.sample() for _ in range(self.mu)]
        if self.warm_start_strategy == WarmStartStrategy.BEST:
            # Add current best individual into starting population
            self.xstart[0] = sorted(self.pop, key=lambda x: x.fitness, reverse=False)[0]
        elif self.warm_start_strategy == WarmStartStrategy.LAST:
            self.xstart[0] = xlast
        self.initialize()
