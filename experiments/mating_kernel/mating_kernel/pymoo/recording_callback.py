from mating_kernel.problems.mo_tracking import MOTrackingMixin
from mating_kernel.pymoo.parser.recording_parser import RecordingParser
from pymoo.core.callback import Callback


def recursive_getattr(obj, attr, default=None):
    for part in attr.split("."):
        obj = getattr(obj, part, default)
        if obj is None:
            return None
    return obj


class RecordingCallback(Callback):
    def __init__(
        self,
        recording_parsers: list[RecordingParser],
        recording_attributes: list[str] | None = None,
    ):
        super().__init__()
        self.recording_parsers = recording_parsers
        self.recording_attributes = (
            recording_attributes if recording_attributes is not None else []
        )
        for parser in recording_parsers:
            self.recording_attributes.extend(parser.record_retrieval)
        self.recording_attributes = list(set(self.recording_attributes))
        if not self.recording_attributes:
            raise ValueError(
                "No recording attributes specified. Please provide a list of recording attributes or at least one RecordingParser."
            )

        self.data["record_stats"] = []

    def notify(self, algorithm):
        records_dict = {}
        for attr in self.recording_attributes:
            operator = recursive_getattr(algorithm, attr)
            if operator is not None and hasattr(operator, "retrieve_records"):
                records = operator.retrieve_records()
                records_dict[attr] = records
        parse_results = {}
        for parser in self.recording_parsers:
            parse_results[parser.__class__.__name__] = parser.parse(records_dict)
        result_dict = parse_results if len(parse_results) > 0 else records_dict
        if isinstance(algorithm.problem.static_problem, MOTrackingMixin):
            result_dict["archive_stats"] = (
                algorithm.problem.static_problem.get_archive_stats()
            )
        self.data["record_stats"].append(result_dict)
