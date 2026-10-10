import numpy as np
import pytest
from typing import Sequence
from pymoo.core.individual import Individual
from pymoo.core.population import Population
from pymoo.core.problem import Problem
from unittest.mock import patch

from mating_kernel.pymoo.pref_model.ornament_pref_model import (
    OrnamentPreferenceModel,
    OrnamentParser,
    Scorer,
)
from mating_kernel.pymoo.mating.preference_mating import PreferenceMating


class DummyOrnamentParser(OrnamentParser):
    def __init__(self, n_dim: int = 1, relative: bool = False):
        assert n_dim > 0, "n_dim must be a positive integer"
        self._n_dim = n_dim
        self._relative = relative

    @property
    def n_dim(self) -> int:
        return self._n_dim

    @property
    def relative(self) -> bool:
        return self._relative

    @property
    def xl(self) -> np.ndarray:
        return np.array([0.0] * self.n_dim)

    @property
    def xu(self) -> np.ndarray:
        return np.array([1.0] * self.n_dim)

    def _ornament(
        self,
        mate: Individual,
        parents: Sequence[Individual] | None = None,
        population: Population | None = None,
        problem: Problem | None = None,
    ) -> np.ndarray:
        # Return a dummy ornament value
        return np.array([0.5] * self.n_dim)


@pytest.mark.parametrize("n_dim", [1, 3])
@pytest.mark.parametrize("relative", [True, False])
def test_ornament_parser(n_dim, relative):
    ind = Individual()
    parents = [Individual(), Individual()]
    population = Population([Individual() for _ in range(5)])
    problem = Problem(n_var=2, n_obj=1)
    # Create a dummy ornament parser
    ornament_parser = DummyOrnamentParser(n_dim=n_dim, relative=relative)

    o_vals = ornament_parser.ornament(
        mate=ind, parents=parents, population=population, problem=problem
    )
    assert isinstance(o_vals, np.ndarray)
    assert o_vals.shape == (n_dim,)
    if relative:
        assert not hasattr(ind, "ornament") or ind.ornament is None
    else:
        assert hasattr(ind, "ornament") and ind.ornament is not None
        assert np.allclose(ind.ornament, o_vals)
    with patch.object(
        ornament_parser, "_ornament", wraps=ornament_parser._ornament
    ) as mock_ornament:
        o_vals2 = ornament_parser.ornament(
            mate=ind, parents=parents, population=population, problem=problem
        )
        assert np.allclose(o_vals, o_vals2)
        if relative:
            # If relative, _ornament should be called again
            mock_ornament.assert_called_once()
        else:
            # If not relative, _ornament should not be called again
            mock_ornament.assert_not_called()


class DummyScorer(Scorer):
    def score(
        self,
        mate: Individual,
        parents: Sequence[Individual] | None = None,
        population: Population | None = None,
        problem: Problem | None = None,
    ) -> float:
        # Return a dummy score based on the mate's ornament values

        return np.random.uniform(0, 1)  # Return a random score for testing


def test_scorer():
    ind = Individual()
    parents = [Individual(), Individual()]
    population = Population([Individual() for _ in range(5)])
    problem = Problem(n_var=2, n_obj=1)
    scorer = DummyScorer()

    score = scorer.score(
        mate=ind, parents=parents, population=population, problem=problem
    )
    assert isinstance(score, float)
    assert 0 <= score <= 1


@pytest.mark.parametrize("n_dim", [1, 3])
@pytest.mark.parametrize("relative", [True, False])
def test_ornament_preference_model(n_dim, relative):
    ornament_parser = DummyOrnamentParser(n_dim=n_dim, relative=relative)
    scorer = DummyScorer()
    model = OrnamentPreferenceModel(ornament_parser=ornament_parser, scorer=scorer)
    assert model.n_preferences == n_dim
    assert model.pref_inheritance is True
    assert model.ornament_parser is ornament_parser
    assert model.scorer is scorer
    assert model.xl.shape == (n_dim,)
    assert model.xu.shape == (n_dim,)

    ind = Individual()
    parents = [Individual(), Individual()]
    population = Population([Individual() for _ in range(5)])
    problem = Problem(n_var=2, n_obj=1)

    # Test initialization of individual
    for m in parents + [ind]:
        model.init_individual(m, random_state=np.random.default_rng())
    assert hasattr(ind, "pref") and ind.pref is not None
    assert ind.pref.shape == (model.n_preferences,)

    # Test evaluation of mate options
    mate_options = [Individual() for _ in range(3)]
    scores = model.evaluate(
        parents=parents,
        mate_options=mate_options,
        problem=problem,
        pop=population,
        random_state=np.random.default_rng(),
    )
    assert isinstance(scores, Sequence)
    assert len(scores) == len(mate_options)

    assert all(isinstance(score, float) for score in scores)


def test_mating_integration():
    ornament_parser = DummyOrnamentParser(n_dim=2, relative=True)
    scorer = DummyScorer()
    model = OrnamentPreferenceModel(ornament_parser=ornament_parser, scorer=scorer)
    assert model.pref_inheritance is True

    from pymoo.operators.selection.tournament import TournamentSelection
    from pymoo.operators.crossover.sbx import SBX
    from pymoo.operators.mutation.pm import PM

    from pymoo.problems import get_problem

    from .. import create_pop, dummy_comp

    mating = PreferenceMating(
        selection=TournamentSelection(func_comp=dummy_comp),
        crossover=SBX(prob=0.9, eta=15),
        mutation=PM(eta=20),
        pref_model=model,
        num_candidates=3,
        candidate_pressure=2,
    )
    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)
    off = mating.do(problem, pop, n_offsprings=5, random_state=np.random.default_rng())

    assert len(off) == 5
    for mate in off:
        assert hasattr(mate, "pref") and mate.pref is not None
        assert mate.pref.shape == (model.n_preferences,)
    for ind in pop:
        assert hasattr(ind, "pref") and ind.pref is not None
        assert ind.pref.shape == (model.n_preferences,)
