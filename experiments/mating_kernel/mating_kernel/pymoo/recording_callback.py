from pymoo.core.callback import Callback

from mating_kernel.problems.mo_tracking import MOTrackingMixin
from mating_kernel.pymoo.parser.recording_parser import RecordingParser, get_record_vars


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
        self.record_retrieval_keys = (
            recording_attributes if recording_attributes is not None else []
        )
        self.record_retrieval_keys.extend(
            get_record_vars(parsers=recording_parsers, attribute="record_retrieval")
        )
        self.record_retrieval_keys = list(set(self.record_retrieval_keys))
        if not self.record_retrieval_keys:
            raise ValueError(
                "No recording attributes specified. Please provide a list of recording attributes or at least one RecordingParser."
            )
        # For convenience, also store the record_args and record_attributes keys from the parsers
        self.record_arg_keys = get_record_vars(
            parsers=recording_parsers, attribute="record_args"
        )
        self.record_attribute_keys = get_record_vars(
            parsers=recording_parsers, attribute="record_attributes"
        )

        self.data["record_stats"] = []

    def notify(self, algorithm):
        records_dict = {}
        for attr in self.record_retrieval_keys:
            operator = recursive_getattr(algorithm, attr)
            if operator is not None and hasattr(operator, "retrieve_records"):
                records = operator.retrieve_records()
                records_dict[attr] = records
        parse_results = {}
        for parser in self.recording_parsers:
            parse_results[parser.__class__.__name__] = parser.parse(records_dict)
        result_dict = parse_results if len(parse_results) > 0 else records_dict
        if hasattr(algorithm.problem, "static_problem") and isinstance(
            algorithm.problem.static_problem, MOTrackingMixin
        ):
            result_dict["archive_stats"] = [
                algorithm.problem.static_problem.get_archive_stats()
            ]  # Adding as a list to maintain consistency with other recorded attributes
        self.data["record_stats"].append(result_dict)
