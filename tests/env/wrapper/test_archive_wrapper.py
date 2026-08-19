from jaix.env.wrapper.archive_wrapper import (
    ArchiveWrapper,
    ArchiveWrapperConfig,
    EnvironmentStepEntry,
)
import gymnasium as gym
from jaix.env.utils.archive.bin_archive import (
    BinArchive,
    BinArchiveConfig,
    BinArchiveEntry,
)
from jaix.env.utils.archive.rv_binning_strategy import (
    RVBinningStrategy,
    RVBinningStrategyConfig,
)
import numpy as np
import pytest
from ..utils.archive.test_archive import DummyArchive, DummyArchiveEntry
from typing import Any
from jaix.env.utils.archive.archive import ArchiveEntry


class DummyEnvironmentStepEntry(
    EnvironmentStepEntry[tuple[np.ndarray, float]], BinArchiveEntry[np.ndarray]
):
    def __init__(
        self,
        action: np.ndarray,
        obs: np.ndarray,
        reward: float,
        terminated: bool,
        truncated: bool,
        info: dict,
    ):
        super().__init__(action, obs, reward, terminated, truncated, info)

    def parse(self) -> tuple[np.ndarray, float]:
        return self.obs, self.reward


class DummyStepArchivEntry(
    DummyArchiveEntry,
    EnvironmentStepEntry[dict[str, Any]],
):
    def __init__(
        self,
        action: Any,
        obs: Any,
        reward: float,
        terminated: bool,
        truncated: bool,
        info: dict,
    ):
        DummyArchiveEntry.__init__(self, sample=obs, fitness=reward)
        EnvironmentStepEntry.__init__(
            self, action, obs, reward, terminated, truncated, info
        )


class DummyStepArchive(DummyArchive):
    @property
    def archive_entry_type(self) -> type[ArchiveEntry]:
        return DummyStepArchivEntry


def do_init(replace_reward=True):
    bin_archive_config = BinArchiveConfig(
        n_bins=5,
        max_fitness=10.0,
        binning_strategy=RVBinningStrategy,
        binning_config=RVBinningStrategyConfig(),
        archive_entry_class=DummyEnvironmentStepEntry,
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
    assert wrapped_env.environment_step_class is DummyEnvironmentStepEntry
    return wrapped_env


def test_init():
    wrapped_env = do_init()
    assert isinstance(wrapped_env.archive, BinArchive)
    assert wrapped_env.replace_reward is True
    assert wrapped_env.passthrough is True


@pytest.mark.parametrize("replace_reward", [True, False])
def test_step(replace_reward):
    wrapped_env = do_init(replace_reward=replace_reward)
    wrapped_env.reset()

    obs, reward, term, trunc, info = wrapped_env.step(wrapped_env.action_space.sample())
    assert isinstance(obs, np.ndarray)
    assert isinstance(reward, float)
    assert isinstance(term, bool)
    assert isinstance(trunc, bool)
    assert isinstance(info, dict)

    assert wrapped_env.archive is not None
    assert wrapped_env.archive.coverage >= 0.0
    assert wrapped_env.archive.coverage <= 1.0


def test_simple_archive():
    env = gym.make("MountainCar-v0", render_mode="rgb_array")
    config = ArchiveWrapperConfig(
        archive_class=DummyStepArchive,
        archive_config=None,
        replace_reward=True,
        passthrough=True,
    )
    wrapped_env = ArchiveWrapper(config, env, max_size=10)
    assert isinstance(wrapped_env.archive, DummyStepArchive)

    wrapped_env.reset()

    obs, reward, term, trunc, info = wrapped_env.step(wrapped_env.action_space.sample())
    assert isinstance(obs, np.ndarray)
    assert isinstance(reward, float)
    assert isinstance(term, bool)
    assert isinstance(trunc, bool)
    assert isinstance(info, dict)

    assert wrapped_env.archive.size == 1
