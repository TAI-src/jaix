from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize

from mating_kernel.problems.mo_tracking import make_tracked
from mating_kernel.pymoo.parser.population_parser import PopulationParser
from mating_kernel.pymoo.parser.reproduction_parser import ReproductionParser
from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper
from mating_kernel.pymoo.recordable_object import make_recordable
from mating_kernel.pymoo.recording_callback import RecordingCallback


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

    parsers = [
        ReproductionParser(ideal=problem.ideal_point),
        PopulationParser(ideal=problem.ideal_point),
    ]
    callback = RecordingCallback(recording_parsers=parsers)

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
    for i, entry in enumerate(records):
        assert "ReproductionParser" in entry
        rep_data = entry["ReproductionParser"]
        if i == 0:  # first generation, no offspring yet
            assert len(rep_data) == 0
        else:
            assert len(rep_data) > 0
            assert "o_F" in rep_data[0]
        assert "PopulationParser" in entry
        pop_data = entry["PopulationParser"]
        assert "n_gen_mean" in pop_data[0]
        assert "archive_stats" in entry
        archive_stats = entry["archive_stats"]
        assert "size" in archive_stats
        assert archive_stats["size"] > 0
        assert archive_stats["fevals"] == (i + 1) * algorithm.pop_size
