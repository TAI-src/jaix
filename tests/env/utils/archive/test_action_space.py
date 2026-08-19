from jaix.env.utils.archive.action_space import (
    ArchiveActionSpace,
    UniformCrossoverActionSpace,
    UniformCrossoverActionSpaceConfig,
)
from ttex.config import ConfigurableObjectFactory as COF
from .test_archive import DummyArchive
import numpy as np
from gymnasium.spaces import MultiDiscrete


def create_dummy_archive_with_samples(
    num_samples: int, max_size: int = 10
) -> DummyArchive:
    archive = DummyArchive(max_size=max_size)
    for i in range(num_samples):
        archive.add(sample=np.array([0, int(i)]), fitness=float(i))
    return archive


def test_archive_action_space_translate():
    archive = create_dummy_archive_with_samples(num_samples=5, max_size=10)
    act_space = MultiDiscrete([archive.max_size] * 3)
    action_space = ArchiveActionSpace(archive, act_space)

    action = [0, 2, 1]
    picked = action_space.translate(action)
    for i, p in enumerate(picked):
        assert p == archive.get(
            action[i]
        ), f"Picked sample {p} does not match expected {archive.get(action[i])}"
        assert (
            p[1] == archive.added_samples[action[i]][1]
        ), f"Picked fitness {p[1]} does not match expected {archive.added_samples[action[i]][1]}"


def test_sample():
    archive = create_dummy_archive_with_samples(num_samples=5, max_size=10)
    act_space = MultiDiscrete([archive.max_size] * 3)
    action_space = ArchiveActionSpace(archive, act_space)
    sampled_action = action_space.sample()

    assert len(sampled_action) == 3, "Sampled action should have length 3"
    assert action_space.action_space.contains(
        sampled_action
    ), "Sampled action should be contained in the action space"


def test_contains():
    archive = create_dummy_archive_with_samples(num_samples=5, max_size=10)
    act_space = MultiDiscrete([archive.max_size] * 3)
    action_space = ArchiveActionSpace(archive, act_space)

    action = [0, 2, 1]
    assert action_space.contains(
        action
    ), "Action should be contained in the action space"

    action = [0, 5, 1]  # only added 5 samples, index 5 is out of bounds
    assert not action_space.contains(
        action
    ), "Action should not be contained in the action space"

    action = [-1, 11, 0]  # truly out of bounds indices
    assert not action_space.contains(
        action
    ), "Action should not be contained in the action space"


def test_uniform_crossover_action_space_translate():
    config = UniformCrossoverActionSpaceConfig(num_parents=3)
    archive = create_dummy_archive_with_samples(num_samples=5, max_size=10)
    action_space = COF.create(UniformCrossoverActionSpace, config, archive=archive)

    action = [0, 2, 1]
    child = action_space.translate(action)

    assert child[0] == 0
    assert child[1] == np.mean(action)


def test_simple():
    config = UniformCrossoverActionSpaceConfig(num_parents=1)
    archive = DummyArchive(max_size=5)
    for i in range(5):
        archive.add(sample=i, fitness=float(i))
    action_space = COF.create(UniformCrossoverActionSpace, config, archive=archive)
    for i in range(5):
        trans_act = action_space.translate([i])
        assert trans_act == i
