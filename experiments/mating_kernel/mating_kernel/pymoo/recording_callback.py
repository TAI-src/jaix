from copy import deepcopy

from pymoo.core.callback import Callback


def recursive_getattr(obj, attr, default=None):
    for part in attr.split("."):
        obj = getattr(obj, part, default)
        if obj is None:
            return None
    return obj


class RecordingCallback(Callback):
    def __init__(self, recording_attributes: list[str] | None = None):
        super().__init__()
        if recording_attributes is None:
            self.recording_attributes = ["problem", "mating.selection"]
        else:
            self.recording_attributes = recording_attributes

        self.data["record_stats"] = []

    def notify(self, algorithm):
        records_dict = {}
        for attr in self.recording_attributes:
            operator = recursive_getattr(algorithm, attr)
            if operator is not None and hasattr(operator, "retrieve_records"):
                records = operator.retrieve_records()
                records_dict[attr] = deepcopy(records)
        self.data["record_stats"].append(records_dict)
