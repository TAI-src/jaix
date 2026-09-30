from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize

from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper
from mating_kernel.pymoo.recording_callback import RecordingCallback

from mating_kernel.problems.mo_tracking import make_tracked
from mating_kernel.pymoo.parser.reproduction_parser import ReproductionParser
from mating_kernel.pymoo.recordable_object import make_recordable


def test_recording_callback():
    # create the algorithm object
    RecordedNSGA2 = make_recordable(NSGA2)
    algorithm = RecordedNSGA2(
        pop_size=5,
        record_args=ReproductionParser.record_args,
        record_attributes=ReproductionParser.record_attributes,
    )

    tracked_REProblem = make_tracked(REProblem)
    problem = tracked_REProblem(REProblemConfig(), inst=0)

    pymoo_problem = PymooProblemWrapper(problem)

    parser = ReproductionParser()
    callback = RecordingCallback(recording_parsers=[parser])

    # execute the optimization
    minimize(
        pymoo_problem,
        algorithm,
        seed=1,
        termination=("n_gen", 5),
        callback=callback,
    )
    records = callback.data["record_stats"]
    assert len(records) == 5
    print(records[1])
