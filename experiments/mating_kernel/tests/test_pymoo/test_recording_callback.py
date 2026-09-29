from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from mating_kernel.pymoo.recording_callback import RecordingCallback
from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper
from .test_do_recorder import RecordedTournamentSelection, dummy_comp


def test_recording_callback():
    # create the algorithm object
    selection = RecordedTournamentSelection(func_comp=dummy_comp)
    algorithm = NSGA2(pop_size=92, selection=selection)

    problem = REProblem(REProblemConfig(), inst=0)
    pymoo_problem = PymooProblemWrapper(problem)

    callback = RecordingCallback()

    # execute the optimization
    minimize(
        pymoo_problem,
        algorithm,
        seed=1,
        termination=("n_gen", 5),
        callback=callback,
    )
    records = callback.data["record_stats"]
    assert len(records) == 5  # 5 generations
    for record in records:
        assert "problem" in record
        assert "mating.selection" in record
        assert len(record["problem"]) == 92  # 92 individuals per generation
    assert (
        len(records[1]["mating.selection"]) > 0
    )  # There should be some records for the selection operator
