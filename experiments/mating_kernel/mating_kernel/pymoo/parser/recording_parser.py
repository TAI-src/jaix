from abc import ABC, abstractmethod
from typing import ClassVar


class RecordingParser(ABC):
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        required_attributes = ["record_args", "record_attributes", "record_retrieval"]
        for attr in required_attributes:
            if not hasattr(cls, attr):
                raise NotImplementedError(
                    f"{cls.__name__} must define the '{attr}' class attribute."
                )

    record_args: ClassVar[list[str]]
    record_attributes: ClassVar[list[str]]
    record_retrieval: ClassVar[list[str]]

    @abstractmethod
    def parse(self, data: dict[str, list]) -> list[dict]: ...
