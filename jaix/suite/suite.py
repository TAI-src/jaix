from abc import abstractmethod
from enum import Enum
from typing import Optional


class AggType(Enum):
    NONE = 0
    INST = 1


class Suite:
    @abstractmethod
    def get_envs(self, agg_type: AggType = AggType.NONE, seed: Optional[int] = None):
        raise NotImplementedError()
