from jaix.env.singular import LJClustEnvironment, LJClustEnvironmentConfig
import pytest
from jaix.env.utils.ase import LJClustAdapterConfig
import numpy as np


@pytest.fixture(scope="session", autouse=True)
def skip_remaining_tests():
    try:
        import ase  # noqa: F401

        assert LJClustAdapterConfig is not None
    except ImportError:
        assert LJClustAdapterConfig is None
        pytest.skip(
            "Skipping LJClust tests. If this is unexpected, check that the ase extra is installed."
        )


@pytest.fixture
def config(skip_remaining_tests):
    return LJClustEnvironmentConfig(
        ljclust_adapter_config=LJClustAdapterConfig(target_dir="./tmp_data"),
        target_accuracy=0.0,
    )


@pytest.fixture
def env(config):
    env = LJClustEnvironment(config, func=0, inst=0)
    return env


def test_init(env):
    assert isinstance(
        env, LJClustEnvironment
    ), "Environment is not an instance of LJClustEnvironment."
    assert env.action_space.shape == (
        env.adapter.num_atoms * 3,
    ), "Action space shape is incorrect."
    assert env.observation_space.shape == (1,), "Observation space shape is incorrect."
    assert env.best_so_far == np.inf, "Best so far should be initialized to infinity."
    assert env.adapter is not None, "Adapter should be initialized."


def test_info(config):
    info = LJClustEnvironment.info(config)
    assert "num_funcs" in info, "Number of functions information is missing."
    assert "num_insts" in info, "Number of instances information is missing."


def test_reset(env):
    env.reset(options={"online": True})
    assert env.best_so_far == np.inf, "Best so far should be reset to infinity."
    info = env._get_info()
    assert "species" in info, "Species information is missing."
    assert "num_atoms" in info, "Number of atoms information is missing."
    assert "box_length" in info, "Box length information is missing."
    assert "min_val" in info, "Minimum value information is missing."
    assert "best_so_far" in info, "Best so far information is missing."


def test_step(env):
    env.reset(options={"online": True})
    for _ in range(10):
        action = env.action_space.sample()
        obs, reward, done, truncated, info = env.step(action)

        assert (
            obs in env.observation_space
        ), "Observation is not in the observation space."
        assert reward is None, "Reward should be None in this environment."
        assert not done, "Environment should not be done after one step."
        assert not truncated, "Environment should not be truncated after one step."
        assert "species" in info, "Species information is missing in step info."

        # Check if the best_so_far is updated
        assert (
            env.best_so_far <= obs
        ), "Best so far should be updated to the minimum reward."


def test_stop(env):
    env.reset(options={"online": True})
    env.best_so_far = env.adapter.min_val + env.target_accuracy + 1e-6
    assert (
        not env.stop()
    ), "Environment should not stop before reaching target accuracy."

    env.best_so_far = env.adapter.min_val + env.target_accuracy
    assert env.stop(), "Environment should stop after reaching target accuracy."


def test_instance_seeding(env):
    # TODO: Implement a test for instance test_instance_seeding
    # This is a placeholder for the actual test.
    pass


def test_render(env):
    # Test rendering functionality
    env.reset(options={"online": True})
    try:
        env.render()
    except NotImplementedError:
        pass  # Rendering is not implemented, which is expected in this case.
    except Exception as e:
        pytest.fail(f"Unexpected error during rendering: {e}")
