from jaix.env.wrapper.archive_wrapper import ArchiveWrapper, ArchiveWrapperConfig
import gymnasium as gym
from jaix.env.utils.archive.bin_archive import BinArchive, BinArchiveConfig
from jaix.env.utils.archive.rv_binning_strategy import (
    RVBinningStrategy,
    RVBinningStrategyConfig,
)
import numpy as np
import pytest


def do_init(replace_reward=True):
    bin_archive_config = BinArchiveConfig(
        n_bins=5,
        max_fitness=10.0,
        binning_strategy=RVBinningStrategy,
        binning_config=RVBinningStrategyConfig(),
        np_bin=1,
        coverage_weight=0.7,
    )
    env = gym.make("MountainCar-v0", render_mode="rgb_array")
    config = ArchiveWrapperConfig(
        archive_class=BinArchive,
        archive_config=bin_archive_config,
        replace_reward=replace_reward,
        passthrough=True,
    )
    ref_dirs = np.array([[1.0], [0.5], [0.0]])
    ideal = np.array([0.0])
    nadir = np.array([10.0])
    wrapped_env = ArchiveWrapper(
        config, env, ref_dirs=ref_dirs, ideal=ideal, nadir=nadir
    )
    return wrapped_env


def test_init():
    wrapped_env = do_init()
    assert isinstance(wrapped_env.archive, BinArchive)
    assert wrapped_env.replace_reward is True
    assert wrapped_env.passthrough is True


@pytest.mark.parametrize("replace_reward", [True, False])
def test_step(replace_reward):
    wrapped_env = do_init(replace_reward=replace_reward)
    ref_dirs = np.array([[1.0], [0.5], [0.0]])
    ideal = np.array([0.0])
    nadir = np.array([10.0])
    wrapped_env.reset(ref_dirs=ref_dirs, ideal=ideal, nadir=nadir)

    obs, reward, term, trunc, info = wrapped_env.step(wrapped_env.action_space.sample())
    assert isinstance(obs, np.ndarray)
    assert isinstance(reward, float)
    assert isinstance(term, bool)
    assert isinstance(trunc, bool)
    assert isinstance(info, dict)

    assert wrapped_env.archive is not None
    assert wrapped_env.archive.coverage >= 0.0
    assert wrapped_env.archive.coverage <= 1.0
