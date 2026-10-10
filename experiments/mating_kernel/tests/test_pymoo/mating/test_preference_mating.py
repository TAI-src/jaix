from mating_kernel.pymoo.mating.preference_mating import PreferenceMating
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.selection.tournament import TournamentSelection
from pymoo.problems import get_problem
from pymoo.algorithms.moo.nsga2 import NSGA2, binary_tournament
from pymoo.optimize import minimize
from ..pref_model.test_pref_model import DummyPreferenceModel


def test_integration_with_nsga2():
    selection = TournamentSelection(func_comp=binary_tournament)
    crossover = SBX(prob=0.9, eta=15)
    mutation = PM(eta=20)
    mating = PreferenceMating(
        selection=selection,
        crossover=crossover,
        mutation=mutation,
        pref_model=DummyPreferenceModel(),
        num_candidates=5,
    )

    problem = get_problem("zdt1")
    algorithm = NSGA2(
        mating=mating,
    )
    res = minimize(
        problem,
        algorithm,
        termination=("n_gen", 5),
        seed=42,
        save_history=True,
    )

    assert res.history is not None
    assert len(res.history) == 5  # 5 generations
