from typing import Any, cast


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
        output = cast(Any, super()).do(*args, **kwargs)
        self.records.append({"args": args, "kwargs": kwargs, "output": output})
        return output

    def retrieve_records(self):
        records = self.records
        self.init_record()
        return records
