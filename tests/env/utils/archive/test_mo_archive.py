from .test_entry_scorer import (
    MOEvalEntry,
    get_problem,
    generate_entries,
    CriticalFrontModes,
)
from jaix.env.utils.archive.mo_archive import MOArchive, MOArchiveConfig, KeepDominated
import pytest
from moocore import pareto_rank, is_nondominated
from jaix.env.utils.archive.entry_scorer import HVContributionScorer
from ttex.config import ConfigurableObjectFactory as COF
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)
from jaix.env.singular.ec_env import ECEnvironment, ECEnvironmentConfig
import numpy as np
from unittest.mock import patch
from itertools import product


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

    for i, entry in enumerate(entries):
        assert entry.rank == ranks[i]


def create_env():

    func = COF.create(REProblem, REProblemConfig(), 1)
    env_config = ECEnvironmentConfig(budget_multiplier=1)
    env = COF.create(ECEnvironment, env_config, func, 0, 1)
    return env


@pytest.mark.parametrize("n_samples", [262_144, None])
def test_archive_init(n_samples):
    config = MOArchiveConfig(
        max_size=10,
        archive_entry_class=MOEvalEntry,
        secondary_criterion_class=HVContributionScorer,
        hv_approx_samples=n_samples,
    )
    archive = MOArchive(config, env=create_env())
    assert archive.size == 0
    assert archive.max_size == 10
    assert isinstance(archive.num_refpoints, int)
    archive.nadir_point = np.array([1.0, 1.0])
    entry = MOEvalEntry(np.array([0.0]), np.array([0.5, 0.5]))
    # Force add
    archive.archived_entries.append(entry)
    score = archive.score
    assert np.isclose(
        score, 0.25, atol=1e-2
    ), f"Score {score} not close to expected 0.25"


def test_hv_caching():
    config = MOArchiveConfig(
        max_size=10,
        archive_entry_class=MOEvalEntry,
        secondary_criterion_class=HVContributionScorer,
        hv_approx_samples=262_144,
    )
    archive = MOArchive(config, env=create_env())
    assert np.isnan(archive._hv)
    entry = MOEvalEntry(np.array([0.0]), np.array([0.5, 0.5]))
    with patch(
        "jaix.env.utils.archive.mo_archive.hv_approx", return_value=-1
    ) as mock_hv:
        archive.nadir_point = np.array([1.0, 1.0])
        archive.archived_entries.append(entry)
        score1 = archive.score
        score2 = archive.score
        assert score1 == score2
        # score called twice but cached value used, so hv_approx should be called only once
        mock_hv.assert_called_once()
    assert not np.isnan(archive._hv)
    # Check that the reset works
    archive.archived_entries = [entry]
    assert np.isnan(archive._hv)


variants = {
    "keep_dominated": [KeepDominated.ALL, KeepDominated.NONE, KeepDominated.STRICT],
    "only_new_entries": [True, False],
    "max_size": [None, 5, 10],
}
# all combinations of the above variants
configurations = [
    dict(zip(variants.keys(), values)) for values in product(*variants.values())
]


@pytest.mark.parametrize("config", configurations)
def test_add_entry_to_archive(config):
    mo_config = MOArchiveConfig(
        archive_entry_class=MOEvalEntry,
        secondary_criterion_class=HVContributionScorer,
        max_size=config["max_size"],
        keep_dominated=config["keep_dominated"],
        only_new_entries=config["only_new_entries"],
    )
    env = create_env()
    archive = MOArchive(mo_config, env=env)
    problem = env.unwrapped.func
    _, gen1 = generate_entries(20, problem, critical_front=CriticalFrontModes.MIDDLE)
    # Force add all entries to the archive
    archive.archived_entries.extend(gen1)

    better_ranks, worse_ranks = generate_entries(
        20, problem, critical_front=CriticalFrontModes.MIDDLE
    )
    gen2 = better_ranks + worse_ranks
    archive.add(gen2)
    if config["max_size"] is not None:
        assert archive.size <= config["max_size"]
    entries = archive.get_all()
    points = np.array([entry.objectives for entry in entries])
    if config["keep_dominated"] == KeepDominated.NONE:
        assert is_nondominated(points, keep_weakly=True, maximise=False).all()
    elif config["keep_dominated"] == KeepDominated.STRICT:
        assert is_nondominated(points, keep_weakly=False, maximise=False).all()
    elif (
        config["keep_dominated"] == KeepDominated.ALL and config["max_size"] is not None
    ):
        assert archive.size == config["max_size"]
    if config["only_new_entries"]:
        assert all(entry in gen2 for entry in entries)
        comparison_entries = gen2
    else:
        comparison_entries = gen1 + gen2
    comparison_points = np.array([entry.objectives for entry in comparison_entries])
    ranks_kept = pareto_rank(points)
    ranks_comparison = pareto_rank(comparison_points)
    max_rank_kept = max(ranks_kept)
    assert max_rank_kept <= max(ranks_comparison)
    num_better_comparison = sum(r < max_rank_kept for r in ranks_comparison)
    num_better_kept = sum(r < max_rank_kept for r in ranks_kept)
    assert num_better_kept == num_better_comparison
