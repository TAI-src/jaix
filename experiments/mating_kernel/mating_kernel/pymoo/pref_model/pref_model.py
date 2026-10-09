from abc import ABC, abstractmethod
from typing import Sequence

import numpy as np
from pymoo.core.problem import Problem
from pymoo.core.individual import Individual
from pymoo.core.population import Population


class PreferenceModel(ABC):
    pref_inheritance: bool  # Whether preference variables are inherited by offspring. If True, offspring will inherit the preference variables of their parents. If False, offspring will be initialized with random preference variables.

    def __init__(self, **kwargs):
        pass

    @property
    @abstractmethod
    def n_preferences(self) -> int:
        """Number of preference variables."""
        ...

    @property
    @abstractmethod
    def xl(self) -> np.ndarray:
        """Lower bounds for preference variables."""
        ...

    @property
    @abstractmethod
    def xu(self) -> np.ndarray:
        """Upper bounds for preference variables."""
        ...

    def initialize(self, random_state: np.random.Generator) -> np.ndarray:
        """Generate an initial preference vector."""
        return random_state.uniform(
            self.xl,
            self.xu,
        )

    def init_individual(
        self, ind: Individual, random_state: np.random.Generator | None = None
    ):
        """Initialize an individual with a random preference vector."""
        if random_state is None:
            rng = np.random.default_rng()
            ind.pref = self.initialize(rng)
        else:
            ind.pref = self.initialize(random_state)

    def validate(self, preferences: np.ndarray) -> bool:
        """Check that a preference vector is valid."""
        return (
            preferences.shape == (self.n_preferences,)
            and np.all(preferences >= self.xl)
            and np.all(preferences <= self.xu)
        )

    @abstractmethod
    def _evaluate(
        self,
        parents: Sequence[Individual],
        mate_options: Sequence[Individual],
        problem: Problem,
        pop: Population,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ) -> Sequence[float]:
        """Evaluates mate options.
        Args:
            parents: The parent individuals.
            mate_options: The candidate mate individuals.
        Returns:
            A preference score for each mate option, where higher scores indicate more preferred mates.
        """
        ...

    def evaluate(
        self,
        parents: Sequence[Individual],
        mate_options: Sequence[Individual],
        problem: Problem,
        pop: Population,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ) -> Sequence[float]:
        """Evaluates mate options and returns preference scores."""
        if self.pref_inheritance:
            # If preferences are inherited
            # make sure every individual has preferences. otherwise initialise.
            for ind in pop:
                if not hasattr(ind, "pref"):
                    self.init_individual(ind, random_state=random_state)
        return self._evaluate(
            parents, mate_options, problem, pop, random_state=random_state, **kwargs
        )

    def dummy_problem(self) -> Problem:
        problem = Problem(
            n_var=self.n_preferences,
            n_obj=1,
            n_constr=0,
            xl=self.xl,
            xu=self.xu,
        )
        return problem
