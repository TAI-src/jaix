from abc import ABC, abstractmethod
from collections import defaultdict
from collections.abc import Sequence
from typing import ClassVar

from mating_kernel.pymoo.do_recorder import RecordingConfig


class RecordingParser(ABC):
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        required_attributes = ["record_args", "record_attributes", "record_retrieval"]
        for attr in required_attributes:
            if not hasattr(cls, attr):
                raise NotImplementedError(
                    f"{cls.__name__} must define the '{attr}' class attribute."
                )

    record_args: ClassVar[dict[str, RecordingConfig]]
    record_attributes: ClassVar[dict[str, RecordingConfig]]
    record_retrieval: ClassVar[list[str]]

    @abstractmethod
    def parse(self, data: dict[str, list]) -> list[dict]: ...


def merge_rec_confs(confs: Sequence[RecordingConfig]) -> RecordingConfig:
    merged_conf = RecordingConfig(
        save_args=list({arg for conf in confs for arg in conf.save_args}),
        cpy_args=list({arg for conf in confs for arg in conf.cpy_args}),
        save_kwargs=list({kw for conf in confs for kw in conf.save_kwargs}),
        cpy_kwargs=list({kw for conf in confs for kw in conf.cpy_kwargs}),
        save_output=any(conf.save_output for conf in confs),
        cpy_output=any(conf.cpy_output for conf in confs),
        meta_fields=list({field for conf in confs for field in conf.meta_fields}),
    )
    return merged_conf


def get_record_vars_config(
    parsers: Sequence[RecordingParser | type[RecordingParser]], attribute: str
) -> dict[str, RecordingConfig]:
    assert attribute in [
        "record_args",
        "record_attributes",
    ], f"Invalid attribute: {attribute}. Must be one of 'record_args', 'record_attributes', 'record_retrieval'."
    configs = [getattr(p, attribute) for p in parsers]
    # First collect all the configurations for each variable
    record_vars: dict[str, list[RecordingConfig]] = defaultdict(list[RecordingConfig])
    for conf in configs:
        if isinstance(conf, dict):
            for key, rec_conf in conf.items():
                record_vars[key].append(rec_conf)
    # Now merge the configurations for each variable
    ret_record_vars: dict[str, RecordingConfig] = {
        key: merge_rec_confs(confs) for key, confs in record_vars.items()
    }
    return ret_record_vars


def get_record_vars_strlist(
    parsers: Sequence[RecordingParser | type[RecordingParser]], attribute: str
) -> list[str]:
    assert attribute in [
        "record_retrieval",
    ], f"Invalid attribute: {attribute}. Must be one of 'record_args', 'record_attributes', 'record_retrieval'."
    record_vars = list({var for p in parsers for var in getattr(p, attribute)})
    return record_vars


def get_record_vars(
    parsers: Sequence[RecordingParser | type[RecordingParser]], attribute: str
) -> dict[str, RecordingConfig] | list[str]:
    if attribute in ["record_args", "record_attributes"]:
        return get_record_vars_config(parsers, attribute)
    elif attribute == "record_retrieval":
        return get_record_vars_strlist(parsers, attribute)
    else:
        raise ValueError(
            f"Invalid attribute: {attribute}. Must be one of 'record_args', 'record_attributes', 'record_retrieval'."
        )
