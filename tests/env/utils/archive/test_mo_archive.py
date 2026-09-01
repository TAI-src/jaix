from .test_entry_scorer import (
    MOEvalEntry,
    get_problem,
    generate_entries,
    CriticalFrontModes,
)
from jaix.env.utils.archive.mo_archive import MOArchive, MOArchiveConfig, KeepDominated
import pytest
from moocore import pareto_rank
from jaix.env.utils.archive.entry_scorer import HVContributionScorer
from ttex.config import ConfigurableObjectFactory as COF
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from jaix.env.singular.ec_env import ECEnvironment, ECEnvironmentConfig
import nupy as np


@pytest.mark.parametrize("keep_dominated", list(KeepDominated))
def test_remove_dominated_entries(keep_dominated):
    problem = get_problem()
    better_ranks, worse_ranks = generate_entries(
        20, problem, critical_front=CriticalFrontModes.MIDDLE
    )
    entries = better_ranks + worse_ranks
    entries = MOArchive.remove_dominated(entries, keep_dominated=keep_dominated)
    assert len(entries) <= len(better_ranks)


@pytest.mark.parametrize("critical_front", list(CriticalFrontModes))
def test_crit_rank(critical_front):
    problem = get_problem()
    _, entries = generate_entries(20, problem, critical_front=CriticalFrontModes.FIRST)
    ranks = pareto_rank([entry.objectives for entry in entries])
    if critical_front == CriticalFrontModes.FIRST:
        crit_rank = 0
    elif critical_front == CriticalFrontModes.LAST:
        crit_rank = max(ranks)
    elif critical_front == CriticalFrontModes.MIDDLE:
        crit_rank = max(ranks) // 2
    n_safe = sum(r < crit_rank for r in ranks)
    max_size = n_safe + 1
    calculated_crit_rank = MOArchive.crit_rank(entries, max_size)
    assert calculated_crit_rank == crit_rank


def create_env():

    func = COF.create(REProblem, REProblemConfig(), 1)
    env_config = ECEnvironmentConfig(budget_multiplier=1)
    env = COF.create(ECEnvironment, env_config, func, 0, 1)
    return env


def test_archive_init():
    config = MOArchiveConfig(
        archive_entry_class=MOEvalEntry, secondary_criterion_class=HVContributionScorer
    )
    archive = MOArchive(config, env=create_env())
    assert archive.size == 0
    assert archive.max_size == 10
    entry = MOEvalEntry(np.array([0.0]), np.array([0.5, 0.5]))
