from typing import Any, cast


class DoRecorderMixin:

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.records = []

    def do(self, *args, **kwargs):
        output = cast(Any, super()).do(*args, **kwargs)
        self.records.append({"args": args, "kwargs": kwargs, "output": output})
        return output

    def retrieve_records(self):
        records = self.records
        self.records = []
        return records
