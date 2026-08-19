from abc import ABC, abstractmethod
from typing import Generic, TypeVar

T = TypeVar("T")


class BinningStrategy(ABC, Generic[T]):
    @abstractmethod
    def get_bin(self, sample: T) -> int:
        """
        Given a sample, return the bin index it belongs to.
        """

    @abstractmethod
    def get_k_nearest_bins(self, bidx: int, k: int) -> list:
        """
        Given a bin index, return the indices of the k nearest bins.
        """
