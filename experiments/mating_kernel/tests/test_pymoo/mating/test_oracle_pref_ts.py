from mating_kernel.pymoo.mating.oracle_pref_ts import OraclePrefTournamentSelection
from .test_random_pref_ts import dummy_comp, create_pop
from pymoo.problems import get_problem


def test_eval_mate():
    # Create a dummy problem and population
    problem = get_problem("zdt1")
    pop = create_pop(problem, size=10)

    # Create an instance of OraclePrefTournamentSelection
    selection = OraclePrefTournamentSelection(func_comp=dummy_comp, num_offspring=5)

    # Select two parents from the population
    parent_idx = [0, 1]  # Indices of the parents in the population
    mate_idx = 2  # Index of the mate in the population

    # Evaluate the mate
    survival_rate = selection.eval_mate(
        [pop[i] for i in parent_idx], pop[mate_idx], problem, pop
    )

    # Check that the survival rate is between 0 and 1
    assert 0 <= survival_rate <= 1
