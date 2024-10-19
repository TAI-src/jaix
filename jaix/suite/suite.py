from abc import abstractmethod
from enum import Enum
from typing import Optional


class AggType(Enum):
    NONE = 0
    INST = 1


class Suite:
    @abstractmethod
    def get_envs(self):
        raise NotImplementedError()

    @abstractmethod
    def get_agg_envs(self, agg_type: AggType, seed: Optional[int] = None):
        raise NotImplementedError()
