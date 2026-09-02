from enum import Enum

import numpy as np
import pytest
from jaix.env.utils.archive.entry_scorer import (
    HVContributionScorer,
    ReferenceVectorDistanceScorer,
)
from jaix.env.utils.archive.mo_archive import MOArchiveEntry
from jaix.env.utils.mo_sizing import get_ref_dirs
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from moocore import pareto_rank


def get_problem(inst: int = 0):
    problem = REProblem(REProblemConfig(), inst)
    return problem


class MOEvalEntry(MOArchiveEntry):
    def __init__(self, x: np.ndarray, y: np.ndarray):
        self.x = x
        self.y = y

    def parse(self) -> np.ndarray:
        return self.y


class CriticalFrontModes(Enum):
    FIRST = "first"
    MIDDLE = "middle"
    LAST = "last"


def generate_entries(
    n_samples: int, problem: REProblem, critical_front: CriticalFrontModes
) -> tuple[list[MOEvalEntry], list[MOEvalEntry]]:
    entries = []
    # Seed so we are sure that the later assertions are correct
    rng = np.random.default_rng(0)
    for _ in range(n_samples):
        x = rng.uniform(problem.lower_bounds, problem.upper_bounds)
        y_raw, _ = problem(x)
        entry = MOEvalEntry(x, np.array(y_raw))
        entries.append(entry)

    # compute fronts
    points = np.array([entry.objectives for entry in entries])
    ranks = pareto_rank(points)
    assert max(ranks) >= 2  # Otherwise the modes don't work.
    # sort entries by rank
    sorted_entries = sorted(zip(entries, ranks), key=lambda er: er[1])
    sorted_entries, ranks = zip(*sorted_entries)
    sorted_entries = list(sorted_entries)
    if critical_front == CriticalFrontModes.FIRST:
        critical_rank = 0
    elif critical_front == CriticalFrontModes.MIDDLE:
        critical_rank = int(max(ranks) // 2)
    elif critical_front == CriticalFrontModes.LAST:
        critical_rank = max(ranks)
    ind_last_safe = sum(rank < critical_rank for rank in ranks)
    assert ind_last_safe > 0 or CriticalFrontModes.FIRST == critical_front
    num_critical = sum(rank == critical_rank for rank in ranks)
    assert num_critical > 0
    ind_last_critical = ind_last_safe + num_critical
    return (
        sorted_entries[:ind_last_safe],
        sorted_entries[ind_last_safe:ind_last_critical],
    )


@pytest.mark.parametrize("mode", list(CriticalFrontModes))
def test_hv_contribution_scorer(mode):
    problem = get_problem(0)
    accepted_entries, critical_entries = generate_entries(20, problem, mode)
    scorer = HVContributionScorer(nadir_point=problem.nadir_point)
    scores = scorer.score(critical_entries, accepted_entries)
    assert len(scores) == len(critical_entries)


@pytest.mark.parametrize("mode", list(CriticalFrontModes))
def test_reference_vector_distance_scorer(mode):
    problem = get_problem(0)
    accepted_entries, critical_entries = generate_entries(20, problem, mode)
    # Create reference directions
    ref_dirs = get_ref_dirs(problem.num_objectives, "original")
    scorer = ReferenceVectorDistanceScorer(
        ref_dirs=ref_dirs,
        ideal_point=problem.ideal_point,
        nadir_point=problem.nadir_point,
    )
    scores = scorer.score(critical_entries, accepted_entries)
    assert len(scores) == len(critical_entries)


def test_hv_order():
    accepted_entries = []
    critical_entries = [
        MOEvalEntry(np.array([0.0, 0.0]), np.array([0.5, 0.5])),
        MOEvalEntry(np.array([0.0, 0.0]), np.array([0.7, 0.45])),
    ]
    scorer = HVContributionScorer(nadir_point=np.array([1.0, 1.0]))
    scores = scorer.score(critical_entries, accepted_entries)
    assert scores[0] > scores[1]


def test_ref_vector_distance_order():
    accepted_entries = []
    critical_entries = [
        MOEvalEntry(np.array([0.0, 0.0]), np.array([0.5, 0.5])),
        MOEvalEntry(np.array([0.0, 0.0]), np.array([0.75, 0.35])),
    ]
    ref_dirs = np.array([[1.0, 0.0], [0.0, 1.0]])
    scorer = ReferenceVectorDistanceScorer(
        ref_dirs=ref_dirs,
        nadir_point=np.array([1.0, 1.0]),
        ideal_point=np.array([0.0, 0.0]),
    )
    scores = scorer.score(critical_entries, accepted_entries)
    assert (
        scores[0] < scores[1]
    )  # The second entry is closer to the reference direction [1, 0] than the first entry.
