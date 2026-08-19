from typing import Any

import numpy as np
from gymnasium.spaces import MultiDiscrete, Space
from ttex.config import Config, ConfigurableObject

from jaix.env.utils.archive.archive import Archive


class ArchiveActionSpace(Space):
    def __init__(self, archive: Archive, action_space: Space):
        self.archive = archive
        self.action_space = action_space

    def translate(self, action: np.ndarray) -> Any:
        assert self.action_space.contains(action), "Action is not in the action space"
        # Pick the solutions from the archive based on the action parameters
        picked = [self.archive.get(i) for i in action]
        return picked

    def sample(self, **kwargs) -> np.ndarray:
        # Sample a random action from the action space
        sample = self.action_space.sample(**kwargs)
        return sample

    @property
    def shape(self) -> tuple[int, ...] | None:
        return self.action_space.shape

    def seed(self, seed: int | None = None) -> int | list[int] | dict[str, int]:
        return self.action_space.seed(seed)

    def contains(self, action: np.ndarray) -> bool:
        if self.action_space.contains(action) is False:
            return False
        # Check if the action is in the action space
        picked = self.translate(action)
        return all(p is not None for p in picked)


class UniformCrossoverActionSpaceConfig(Config):
    def __init__(self, num_parents: int = 2):
        Config.__init__(self)
        self.num_parents = num_parents


class UniformCrossoverActionSpace(ArchiveActionSpace, ConfigurableObject):
    config_class = UniformCrossoverActionSpaceConfig

    def __init__(self, config: UniformCrossoverActionSpaceConfig, archive: Archive):
        ConfigurableObject.__init__(self, config)
        assert archive.max_size is not None, "Archive must have a max_size defined"
        ArchiveActionSpace.__init__(
            self, archive, MultiDiscrete([archive.max_size] * self.num_parents)
        )

    def translate(self, action: np.ndarray) -> Any:
        archive_content = super().translate(action)
        parents = []
        for p in archive_content:
            # parents are tuples of (sample, fitness), we only want the sample
            assert p is not None, "Parent is None, cannot perform crossover"
            assert len(p) == 2, "Parent should be a tuple of (sample, fitness)"
            assert self.action_space.contains(
                action
            ), "Action is not in the action space"
            parents.append(p[0])
        offspring = np.mean(parents, axis=0)
        if self.action_space.dtype is not None:
            # try to convert to dtype of underlying action space
            offspring = np.asarray(offspring, dtype=self.action_space.dtype)
        return offspring
