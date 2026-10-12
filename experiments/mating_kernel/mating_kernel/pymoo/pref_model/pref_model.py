from abc import ABC, abstractmethod
from collections.abc import Sequence

import numpy as np
from pymoo.core.individual import Individual
from pymoo.core.population import Population
from pymoo.core.problem import Problem


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
        random_vector = random_state.uniform(self.xl, self.xu, size=self.n_preferences)
        return random_vector

    def init_individual(
        self, ind: Individual, random_state: np.random.Generator | None = None
    ):
        """Initialize an individual with a random preference vector."""
        if not self.pref_inheritance or (hasattr(ind, "pref") and ind.pref is not None):
            return  # Individual already has a preference vector, do not reinitialize
        if random_state is None:
            rng = np.random.default_rng()
            ind.pref = self.initialize(rng)
        else:
            ind.pref = self.initialize(random_state)

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
        scores = self._evaluate(
            parents, mate_options, problem, pop, random_state=random_state, **kwargs
        )
        return scores

    @property
    def dummy_problem(self) -> Problem:
        problem = Problem(
            n_var=self.n_preferences,
            n_obj=1,
            n_constr=0,
            xl=self.xl,
            xu=self.xu,
        )
        return problem
