from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, cast


@dataclass
class RecordingConfig:
    save_args: list[int] = field(default_factory=list)
    cpy_args: list[int] = field(default_factory=list)
    save_kwargs: list[str] = field(default_factory=list)
    cpy_kwargs: list[str] = field(default_factory=list)
    save_output: bool = False
    cpy_output: bool = False
    meta_fields: list[str] = field(default_factory=list)


class DoRecorderMixin:
    def __init__(
        self,
        recording_config: RecordingConfig | None = None,
    ):
        self.recording_config = recording_config or RecordingConfig()

    @property
    def records(self):
        if not hasattr(self, "_records"):
            self.init_record()
        return self._records

    def init_record(self):
        self._records = []
        meta = {
            field: getattr(self, field) for field in self.recording_config.meta_fields
        }
        self._records.append(meta)

    def do(self, *args, **kwargs):
        args_cpy = [
            deepcopy(args[i]) for i in self.recording_config.cpy_args if i < len(args)
        ]
        args_save = [args[i] for i in self.recording_config.save_args if i < len(args)]
        kwargs_cpy = {
            k: deepcopy(v)
            for k, v in kwargs.items()
            if k in self.recording_config.cpy_kwargs
        }
        kwargs_save = {
            k: v for k, v in kwargs.items() if k in self.recording_config.save_kwargs
        }
        output = cast(Any, super()).do(*args, **kwargs)
        output_cpy = deepcopy(output) if self.recording_config.cpy_output else None
        output_save = output if self.recording_config.save_output else None
        self.records.append(
            {
                "args": args_save,
                "args_cpy": args_cpy,
                "kwargs": kwargs_save,
                "kwargs_cpy": kwargs_cpy,
                "output": output_save,
                "output_cpy": output_cpy,
            }
        )
        return output

    def retrieve_records(self):
        records = self.records
        self.init_record()
        return records
