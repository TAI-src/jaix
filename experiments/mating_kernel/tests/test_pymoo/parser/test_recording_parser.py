from typing import ClassVar

import pytest

from mating_kernel.pymoo.parser.recording_parser import RecordingParser


def test_valid():
    # valid test
    class ValidParser(RecordingParser):
        record_args: ClassVar[list[str]] = ["test"]
        record_attributes: ClassVar[list[str]] = []
        record_retrieval: ClassVar[list[str]] = ["a", "b", "c"]

        def parse(self, data: dict[str, list]) -> list[dict]:
            return []

    ValidParser()


@pytest.mark.parametrize(
    "missing_attribute",
    ["record_args", "record_attributes", "record_retrieval"],
)
def test_recording_parser_subclass_requires_each_attribute(missing_attribute):
    attributes = {
        "record_args": ["test"],
        "record_attributes": [],
        "record_retrieval": ["a", "b", "c"],
    }
    del attributes[missing_attribute]

    with pytest.raises(
        NotImplementedError,
        match=f"TestParser must define the '{missing_attribute}' class attribute.",
    ):
        type("TestParser", (RecordingParser,), attributes)
