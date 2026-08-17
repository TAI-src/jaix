from abc import ABC, abstractmethod
from ttex.config import Config, ConfigurableObject
from jaix.env.archive.archive import Archive
import numpy as np
from typing import Any


class ArchiveAction(ABC):
    """
    A class that represents actions on archives, such as a recombination of two solutions in the archive, or a mutation of a solution in the archive.
    """

    @staticmethod
    def pick(archive: Archive, action: np.ndarray) -> Any:
        # Pick the solutions from the archive based on the action parameters
        picked = [archive.get(i) for i in action]
        return picked

    @abstractmethod
    def execute(self, archive: Archive, action: np.ndarray) -> Any:
        """
        Execute the action
        """

    @abstractmethod
    def random_action(self, archive: Archive) -> Any:
        """
        Generate a random action
        """


class UniformCrossoverConfig(Config):
    def __init__(self, num_parents: int = 2):
        Config.__init__(self)
        self.num_parents = num_parents


class UniformCrossover(ArchiveAction, ConfigurableObject):
    config_class = UniformCrossoverConfig

    def execute(self, archive: Archive, action: np.ndarray) -> Any:
        assert len(action) == self.num_parents
        parents = ArchiveAction.pick(archive, action)
        offspring = np.mean(parents, axis=1)
        return offspring

    def random_action(self, archive: Archive) -> Any:
        action = np.random.randint(low=0, high=archive.size, size=self.num_parents)
        return action
