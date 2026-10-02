from typing import Any, cast
from copy import deepcopy


class DoRecorderMixin:
    @property
    def records(self):
        if not hasattr(self, "_records"):
            self.init_record()
        return self._records

    def init_record(self):
        self._records = []
        self._records.append(self.__dict__)

    def do(self, *args, **kwargs):
        args_cpy = deepcopy(args)
        kwargs_cpy = deepcopy(kwargs)
        output = cast(Any, super()).do(*args, **kwargs)
        self.records.append(
            {
                "args": args,
                "args_cpy": args_cpy,
                "kwargs": kwargs,
                "kwargs_cpy": kwargs_cpy,
                "output": output,
                "output_cpy": deepcopy(output),
            }
        )
        return output

    def retrieve_records(self):
        records = self.records
        self.init_record()
        return records
