import gymnasium as gym
from jaix.env.utils.archive.action_space import (
    UniformCrossoverActionSpace,
    UniformCrossoverActionSpaceConfig,
)
from jaix.env.wrapper.archive_action_wrapper import (
    ArchiveActionWrapper,
    ArchiveActionWrapperConfig,
)
from jaix.env.wrapper.archive_wrapper import ArchiveWrapperConfig

from .test_archive_wrapper import DummyStepArchive, DummyStepArchivEntry


def create_dummy_archive_action_wrapper(
    max_size: int = 10, num_samples: int = 5
) -> ArchiveActionWrapper:

    # Create a dummy archive with samples
    archive_config = ArchiveWrapperConfig(
        archive_class=DummyStepArchive,
        archive_config=None,
        replace_reward=True,
        passthrough=True,
    )

    action_space_config = UniformCrossoverActionSpaceConfig(
        crossover_attribute="action", num_parents=1
    )

    wrapper_config = ArchiveActionWrapperConfig(
        archive_wrapper_config=archive_config,
        action_space_class=UniformCrossoverActionSpace,
        action_space_config=action_space_config,
    )
    env = gym.make("Taxi-v4")  # Dummy environment for testing purposes

    # Create the ArchiveActionWrapper
    wrapper = ArchiveActionWrapper(
        wrapper_config, env=env, max_size=max_size
    )  # env is None for testing purposes

    # Add samples to the archive
    for i in range(num_samples):
        sample = env.action_space.sample()  # Sample from the environment's action space
        wrapper.archive.add(
            [
                DummyStepArchivEntry(
                    action=sample,
                    obs=None,
                    reward=float(i),
                    terminated=False,
                    truncated=False,
                    info={},
                )
            ]
        )

    return wrapper


def test_step():
    wrapped_env = create_dummy_archive_action_wrapper(max_size=10, num_samples=12)

    obs, info = wrapped_env.reset()
    n_points_old = wrapped_env.archive.get_archive_stats()["num_points"]
    action = wrapped_env.action_space.sample()
    expected_offspring = wrapped_env.archive.get(action[0])
    expected_attribute = getattr(
        expected_offspring, wrapped_env.action_space.crossover_attribute
    )
    obs, reward, term, trunc, info = wrapped_env.step(action)
    n_points_new = wrapped_env.archive.get_archive_stats()["num_points"]
    assert (
        n_points_new == n_points_old + 1
    ), "Number of points in the archive should have increased by 1"
    assert (
        info["env_action"] == expected_attribute
    ), "The offspring in the info should match the expected offspring"
    entry = wrapped_env.archive.get(9)
    assert (
        entry.action == expected_attribute
    ), "The added sample should match the expected offspring"
    assert (
        entry.reward == reward
    ), "The added fitness should match the reward returned by the step function"
