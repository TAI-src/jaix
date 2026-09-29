from pymoo.core.callback import Callback


class RecordingCallback(Callback):
    def __init__(self, recording_attributes: list[str] | None = None):
        super().__init__()
        if recording_attributes is None:
            self.recording_attributes = ["problem", "selection"]
        else:
            self.recording_attributes = recording_attributes

        self.data["record_stats"] = {attr: [] for attr in self.recording_attributes}

    def notify(self, algorithm):
        for attr in self.recording_attributes:
            operator = getattr(algorithm, attr, None)
            if operator is not None and hasattr(operator, "retrieve_records"):
                records = operator.retrieve_records()
                self.data["record_stats"][attr].extend(records)
