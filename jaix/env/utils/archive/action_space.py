from typing import Any

import numpy as np
from gymnasium.spaces import MultiDiscrete, Space
from ttex.config import Config, ConfigurableObject

from jaix.env.utils.archive.archive import Archive, ArchiveEntry
from ABC import abstractmethod


class ArchiveActionSpace(Space):
    def __init__(self, archive: Archive, action_space: Space):
        self.archive = archive
        self.action_space = action_space

    @abstractmethod
    def translate(self, action: np.ndarray) -> Any:
        """
        Translate the action into the archive's action space. This method should be implemented by subclasses.
        """
        ...

    def pick(self, action: np.ndarray) -> list[ArchiveEntry | None]:
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
        picked = self.pick(action)
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
        assert self.action_space.contains(action), "Action is not in the action space"
        archive_content = self.pick(action)
        parents = [getattr(p, "x", None) for p in archive_content if p is not None]
        assert (
            len(parents) == self.num_parents
        ), "Not enough parents found in the archive"
        assert all(
            isinstance(p, np.ndarray) for p in parents
        ), "All parents must be numpy arrays"

        # Create narrowed list of parents to avoid issues with None values
        parents = [p for p in parents if isinstance(p, np.ndarray)]
        offspring = np.mean(parents, axis=0)
        if self.action_space.dtype is not None:
            # try to convert to dtype of underlying action space
            offspring = np.asarray(offspring, dtype=self.action_space.dtype)
        return offspring
