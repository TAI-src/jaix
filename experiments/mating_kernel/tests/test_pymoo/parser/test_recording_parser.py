from typing import ClassVar

import pytest

from mating_kernel.pymoo.parser.recording_parser import RecordingParser, get_record_vars


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


class DummyParser(RecordingParser):
    record_args: ClassVar[list[str]] = ["arg1", "arg2"]
    record_attributes: ClassVar[list[str]] = ["attr1", "attr2"]
    record_retrieval: ClassVar[list[str]] = ["ret1", "ret2"]

    def parse(self, data: dict[str, list]) -> list[dict]:
        return []


class DummyParser2(RecordingParser):
    record_args: ClassVar[list[str]] = ["arg1", "arg3"]
    record_attributes: ClassVar[list[str]] = ["attr3"]
    record_retrieval: ClassVar[list[str]] = ["ret3"]

    def parse(self, data: dict[str, list]) -> list[dict]:
        return []


attributes = ["record_args", "record_attributes", "record_retrieval"]
parsers = [
    [DummyParser, DummyParser2],
    [DummyParser(), DummyParser2()],
    [DummyParser, DummyParser2()],
]
combs = [(attr, p) for attr in attributes for p in parsers]


@pytest.mark.parametrize("attribute, parsers", combs)
def test_get_record_vars(attribute, parsers):
    expected = {
        "record_args": ["arg1", "arg2", "arg3"],
        "record_attributes": ["attr1", "attr2", "attr3"],
        "record_retrieval": ["ret1", "ret2", "ret3"],
    }
    result = get_record_vars(parsers, attribute)
    assert set(result) == set(expected[attribute])
