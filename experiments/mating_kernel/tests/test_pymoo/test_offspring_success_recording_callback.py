import pandas as pd
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize

from mating_kernel.problems.mo_tracking import make_tracked
from mating_kernel.pymoo.offspring_success_recording_callback import (
    OffspringSuccessRecordingCallback,
)
from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper
from mating_kernel.pymoo.recordable_object import make_recordable


def run_with_callback(pop_size=5, n_gen=6):

    tracked_REProblem = make_tracked(REProblem)
    problem = tracked_REProblem(REProblemConfig(), inst=0)

    callback = OffspringSuccessRecordingCallback(problem)

    # create the algorithm object
    RecordedNSGA2 = make_recordable(NSGA2)
    algorithm = RecordedNSGA2(
        pop_size=pop_size,
        record_args=callback.record_arg_keys,
        record_attributes=callback.record_attribute_keys,
    )

    # execute the optimization
    minimize(
        PymooProblemWrapper(problem),
        algorithm,
        seed=1,
        termination=("n_gen", n_gen),
        callback=callback,
    )
    return callback


def test_offspring_success_recording_callback(pop_size=5, n_gen=6):
    callback = run_with_callback()
    records = callback.data["record_stats"]
    assert len(records) == n_gen


def test_get_stat_gen(pop_size=5, n_gen=6):
    callback = run_with_callback()
    for gen in range(n_gen):
        for stat_name in ["PopulationParser", "archive_stats"]:
            stat = callback.get_stat_gen(gen, stat_name)
            assert isinstance(stat, dict)


def test_get_stat_gen_out_of_bounds(pop_size=5, n_gen=6):
    callback = run_with_callback()
    try:
        callback.get_stat_gen(-1, "offspring_success")
    except IndexError as e:
        assert str(e) == "Generation -1 is out of bounds."
    try:
        callback.get_stat_gen(n_gen, "offspring_success")
    except IndexError as e:
        assert str(e) == f"Generation {n_gen} is out of bounds."


def test_record_stats_property(pop_size=5, n_gen=6):
    callback = run_with_callback()
    record_stats = callback._get_full_record_dicts()
    assert isinstance(record_stats, list)
    for i, record in enumerate(record_stats):
        assert "survived" in record
        assert "o_X" in record
        assert "b_crowding_mean" in record
        assert "a_crowding_mean" in record
        assert "b_size" in record
        assert "a_size" in record
        assert record["a_fevals"] == record["b_fevals"] + pop_size
    assert record_stats[0]["b_fevals"] == pop_size


def test_dirty_flag(pop_size=5, n_gen=6):
    callback = run_with_callback()
    assert callback._dirty is True
    # After retrieving the full record dicts, the dirty flag should still be True
    _ = callback._get_full_record_dicts()
    assert callback._dirty is True
    _ = callback.record_stats
    assert callback._dirty is False


def test_record_stats(pop_size=5, n_gen=6):
    callback = run_with_callback()
    record_stats = callback.record_stats
    assert isinstance(record_stats, pd.DataFrame)
    assert "survived" in record_stats.columns
