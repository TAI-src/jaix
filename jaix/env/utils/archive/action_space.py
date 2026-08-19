import numbers
from abc import ABC, abstractmethod
from typing import Any

import numpy as np
from gymnasium.spaces import MultiDiscrete, Space
from ttex.config import Config, ConfigurableObject

from jaix.env.utils.archive.archive import Archive, ArchiveEntry


class ArchiveActionSpace(Space, ABC):
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

    def sample(self, mask: Any | None = None, probability: Any | None = None) -> Any:
        # Sample a random action from the action space
        sample = self.action_space.sample(mask=mask, probability=probability)
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
    def __init__(self, crossover_attribute: str, num_parents: int = 2):
        Config.__init__(self)
        self.crossover_attribute = crossover_attribute
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
        assert self.action_space.dtype is not None, "Action space dtype is not defined"
        archive_content = self.pick(action)
        parents = [
            getattr(p, self.crossover_attribute)
            for p in archive_content
            if p is not None
        ]
        assert (
            len(parents) == self.num_parents
        ), "Not enough parents found in the archive"

        # Treat numerical values and numpy arrays differently
        if all(isinstance(p, numbers.Number) for p in parents):
            offspring = self.action_space.dtype.type(np.mean(parents))
        elif all(isinstance(p, np.ndarray) for p in parents):
            offspring = np.asarray(
                np.mean(parents, axis=0), dtype=self.action_space.dtype
            )
        else:
            raise TypeError(
                f"Parents must be all numbers or all numpy arrays, got types: {[type(p) for p in parents]}"
            )
        return offspring
