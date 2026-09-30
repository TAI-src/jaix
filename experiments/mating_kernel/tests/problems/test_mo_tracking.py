from jaix.env.utils.archive.mo_archive import KeepDominated
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from mating_kernel.problems.mo_tracking import MOTrackingMixin, make_tracked
from jaix.env.utils.problem.static_problem import StaticProblem
import numpy as np
import pytest


class TrackedREProblem(MOTrackingMixin, REProblem):
    pass


def test_init():
    problem = TrackedREProblem(REProblemConfig(), inst=0)
    assert isinstance(problem, REProblem)
    assert isinstance(problem, MOTrackingMixin)
    assert isinstance(problem, StaticProblem)

    # Check if the archive is created correctly
    assert problem.archive is not None
    assert problem.archive.config.max_size is None
    assert problem.archive.config.keep_dominated == KeepDominated.NONE
    assert problem.archive.config.only_new_entries is False
    assert (
        problem.archive.config.secondary_criterion_class.__name__
        == "ReferenceVectorDistanceScorer"
    )


def test_eval():
    problem = TrackedREProblem(REProblemConfig(), inst=0)

    for i in range(5):
        x = np.random.uniform(problem.lower_bounds, problem.upper_bounds)
        problem(x)
        archive_entries = problem.archive.queued_entries
        assert len(archive_entries) == i + 1


problem_classes = [
    TrackedREProblem,
    make_tracked(REProblem),
]


@pytest.mark.parametrize("problem_class", problem_classes)
def test_get_stats(problem_class):
    problem = problem_class(REProblemConfig(), inst=0)

    # prefill the archive with random samples
    for _ in range(5):
        x = np.random.uniform(problem.lower_bounds, problem.upper_bounds)
        problem(x)

    stats = problem.get_archive_stats()
    assert stats["fevals"] == 5
    assert stats["size"] > 0
