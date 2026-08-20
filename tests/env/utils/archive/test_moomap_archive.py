from ..problem import CobiProblem, CobiProblemConfig
import pytest
from jaix.env.utils.archive.moomap_archive import (
    MoomapArchive,
    MoomapArchiveEntry,
    MoomapArchiveConfig,
)
from jaix.env.singular.ec_env import ECEnvironment, ECEnvironmentConfig
from ttex.config import ConfigurableObjectFactory as COF
import numpy as np


@pytest.fixture(scope="session", autouse=True)
def skip_remaining_tests():
    try:
        import cobi  # noqa: F401

        assert CobiProblem is not None
    except ImportError:
        assert CobiProblem is None
        pytest.skip(
            "Skipping CobiProblem tests. If this is unexpected, check that the jaix_cobi docker image is used."
        )


def create_env():
    config = CobiProblemConfig(
        n_var=3, n_constraints={"Linear": 0, "Quadratic": 0, "Multi": 0}, domain=(-4, 4)
    )
    func = COF.create(CobiProblem, config, 1)
    env_config = ECEnvironmentConfig(budget_multiplier=1)
    env = COF.create(ECEnvironment, env_config, func, 0, 1)
    return env


def test_moomap_archive_entry():
    env = create_env()

    action = env.action_space.sample()
    obs, reward, terminated, truncated, info = env.step(action)
    entry = MoomapArchiveEntry(action, obs, reward, terminated, truncated, info)

    sample, fitness = entry.parse()
    assert np.array_equal(sample, obs)
    assert np.isnan(fitness)  # fitness is not set yet

    ideal_point = env.unwrapped.func.ideal_point
    entry.set_fitness(ideal_point)
    assert not np.isnan(entry.fitness)  # fitness is now set


def test_moomap_archive_init():
    env = create_env()
    archive_config = MoomapArchiveConfig(np_bin=1, coverage_weight=0.5)
    archive = MoomapArchive(archive_config, env)

    assert isinstance(archive, MoomapArchive)
    assert hasattr(archive, "map") and isinstance(archive.map, dict)

    assert archive.coverage == 0.0, "Initial coverage should be 0.0"
    assert np.isnan(archive.fitness), "Initial fitness should be NaN"
    assert np.isnan(archive.score), "Initial score should be NaN"
    archive_stats = archive.get_archive_stats()
    assert (
        archive.num_refpoints == 100
    ), "Number of reference points should be 100 for 2d"
    assert archive_stats["coverage_100"] == 0.0, "Stats coverage should be 0.0"


def test_add_entry_to_moomap_archive():
    env = create_env()
    archive_config = MoomapArchiveConfig(np_bin=1, coverage_weight=0.5)
    archive = MoomapArchive(archive_config, env)

    action = env.action_space.sample()
    obs, reward, terminated, truncated, info = env.step(action)
    entry = MoomapArchiveEntry(action, obs, reward, terminated, truncated, info)
    # Add entry to the archive
    score = archive.add(entry)
    assert np.isnan(score)  # initial add returns NaN because it is about improvement

    # Check if the entry was added
    assert len(archive.map) > 0, "Entry should be added to the archive"
    assert archive.covered_bins == 1
    assert sum(archive.add_counter) == 1
    assert (
        archive.coverage > 0.0
    ), "Coverage should be greater than 0 after adding an entry"


def test_get():
    env = create_env()
    archive_config = MoomapArchiveConfig(np_bin=1, coverage_weight=0.5)
    archive = MoomapArchive(archive_config, env)

    action = env.action_space.sample()
    obs, reward, terminated, truncated, info = env.step(action)
    entry = MoomapArchiveEntry(action, obs, reward, terminated, truncated, info)
    # Add entry to the archive
    archive.add(entry)

    # Get the entry from the archive
    retrieved_entry = archive.get(0)
    assert retrieved_entry is not None, "Retrieved entry should not be None"
