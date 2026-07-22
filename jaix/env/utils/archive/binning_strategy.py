from abc import ABC, abstractmethod
from typing import Any


class BinningStrategy(ABC):
    @abstractmethod
    def get_bin(self, sample: Any) -> int:
        """
        Given a sample, return the bin index it belongs to.
        """
        pass

    @abstractmethod
    def get_k_nearest_bins(self, bidx: int, k: int) -> list:
        """
        Given a bin index, return the indices of the k nearest bins.
        """
        pass
