from jaix.env.singular import LJClustEnvironment, LJClustEnvironmentConfig
import pytest
from jaix.env.utils.ase import LJClustAdapterConfig


@pytest.fixture
def env():
    config = LJClustEnvironmentConfig(
        ljclust_adapter_config=LJClustAdapterConfig(target_dir="./tmp_data"),
        target_accuracy=0.0,
    )
    env = LJClustEnvironment(config, func=0, inst=0)
    return env


def test_init(env):
    assert isinstance(
        env, LJClustEnvironment
    ), "Environment is not an instance of LJClustEnvironment."
    assert env.action_space.shape == (
        env.adapter.num_atoms,
        3,
    ), "Action space shape is incorrect."
    assert env.observation_space.shape == (1,), "Observation space shape is incorrect."
    assert env.best_so_far == np.inf, "Best so far should be initialized to infinity."
    assert env.adapter is not None, "Adapter should be initialized."


def test_reset(env):
    env.reset()
    assert env.best_so_far == np.inf, "Best so far should be reset to infinity."
    info = env._get_info()
    assert "species" in info, "Species information is missing."
    assert "num_atoms" in info, "Number of atoms information is missing."
    assert "box_length" in info, "Box length information is missing."
    assert "min_val" in info, "Minimum value information is missing."
    assert "best_so_far" in info, "Best so far information is missing."


def test_step(env):
    env.reset()
    action = env.action_space.sample()
    obs, reward, done, truncated, info = env.step(action)

    assert obs in env.observation_space, "Observation is not in the observation space."
    assert isinstance(reward, float), "Reward should be a float."
    assert not done, "Environment should not be done after one step."
    assert not truncated, "Environment should not be truncated after one step."

    # Check if the best_so_far is updated
    assert (
        env.best_so_far >= obs
    ), "Best so far should be updated to the minimum reward."
