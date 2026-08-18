from typing import Any

import numpy as np
from gymnasium.spaces import MultiDiscrete
from ttex.config import Config, ConfigurableObject

from jaix.env.utils.archive.archive import Archive


class ArchiveActionSpace(MultiDiscrete):
    def __init__(self, archive: Archive, num_choices: int, **kwargs):
        self.archive = archive
        self.num_choices = num_choices
        assert self.archive.max_size is not None, "Archive must have a max_size defined"
        super().__init__([self.archive.max_size] * self.num_choices, **kwargs)

    def translate(self, action: np.ndarray) -> Any:
        assert super().contains(action), "Action is not in the action space"
        # Pick the solutions from the archive based on the action parameters
        picked = [self.archive.get(i) for i in action]
        return picked

    def sample(
        self, mask: Any | None = None, probability: Any | None = None
    ) -> np.ndarray:
        # Sample a random action from the action space
        sample = super().sample(mask, probability)
        return self.translate(sample)

    def contains(self, action: np.ndarray) -> bool:
        if super().contains(action) is False:
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

    def __init__(
        self, config: UniformCrossoverActionSpaceConfig, archive: Archive, **kwargs
    ):
        ConfigurableObject.__init__(self, config)
        super().__init__(archive, self.num_parents, **kwargs)

    def translate(self, action: np.ndarray) -> Any:
        archive_content = super().translate(action)
        parents = []
        for p in archive_content:
            # parents are tuples of (sample, fitness), we only want the sample
            assert p is not None, "Parent is None, cannot perform crossover"
            assert len(p) == 2, "Parent should be a tuple of (sample, fitness)"
            assert isinstance(p[0], np.ndarray), "Parent sample should be a numpy array"
            parents.append(p[0])
        offspring = np.mean(parents, axis=0)
        return offspring
