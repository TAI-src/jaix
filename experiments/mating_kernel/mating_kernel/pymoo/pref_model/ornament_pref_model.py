from abc import ABC, abstractmethod
import numpy as np
from pymoo.core.individual import Individual
from pymoo.core.population import Population
from pymoo.core.problem import Problem

from typing import Sequence
from mating_kernel.pymoo.pref_model.pref_model import PreferenceModel


class OrnamentParser(ABC):
    @property
    @abstractmethod
    def n_dim(self) -> int:
        """Number of dimensions in the ornament space."""
        ...

    @property
    @abstractmethod
    def relative(self) -> bool:
        """Whether the ornament values are relative to parents or the population. If True, the ornament values are computed relative to the parents or the population. If False, the ornament values are computed independently of the parents or the population."""
        ...

    @property
    @abstractmethod
    def xl(self) -> np.ndarray: ...

    @property
    @abstractmethod
    def xu(self) -> np.ndarray: ...

    @abstractmethod
    def _ornament(
        self,
        mate: Individual,
        parents: Sequence[Individual] | None = None,
        population: Population | None = None,
        problem: Problem | None = None,
    ) -> np.ndarray:
        """Extract ornament values for the mate
        potentially relative to parents and/or the whole population
        """
        ...

    def ornament(
        self,
        mate: Individual,
        parents: Sequence[Individual] | None = None,
        population: Population | None = None,
        problem: Problem | None = None,
    ) -> np.ndarray:
        """Extract ornament values for the mate
        potentially relative to parents and/or the whole population
        """
        if hasattr(mate, "ornament") and mate.ornament is not None:
            ornament_values = mate.ornament
        elif self.relative:
            ornament_values = self._ornament(
                mate=mate, parents=parents, population=population, problem=problem
            )
        else:
            ornament_values = self._ornament(
                mate=mate, parents=None, population=None, problem=problem
            )
            setattr(
                mate, "ornament", ornament_values
            )  # Store ornament values in mate for later use since they are not relative to others
        return ornament_values


class Scorer(ABC):
    @abstractmethod
    def score(
        self, ornament_values: np.ndarray, parent_preferences: Sequence[np.ndarray]
    ) -> float:
        """Score a mate based on its ornament values
        and the preferences of the previously chosen parents
        """
        ...


class OrnamentPreferenceModel(PreferenceModel):
    pref_inheritance = True  # Preferences are inherited from parents

    def __init__(self, ornament_parser: OrnamentParser, scorer: Scorer, **kwargs):
        super().__init__(**kwargs)
        self.ornament_parser = ornament_parser
        self.scorer = scorer

    @property
    def n_preferences(self) -> int:
        return self.ornament_parser.n_dim

    @property
    def xl(self) -> np.ndarray:
        return self.ornament_parser.xl

    @property
    def xu(self) -> np.ndarray:
        return self.ornament_parser.xu

    def _evaluate(
        self,
        parents: Sequence[Individual],
        mate_options: Sequence[Individual],
        problem: Problem,
        pop: Population,
        random_state: np.random.Generator | None = None,
        **kwargs,
    ) -> Sequence[float]:
        scores = []
        parent_preferences = [parent.pref for parent in parents]
        for mate in mate_options:
            ornament_values = self.ornament_parser.ornament(
                mate=mate, parents=parents, population=pop, problem=problem
            )
            score = self.scorer.score(ornament_values, parent_preferences)
            scores.append(score)

        return scores
