from typing import Any, Dict, Optional
from uuid import uuid4
from dataclasses import dataclass


@dataclass
class Artifact:
    name: str
    local_path: str
    artifact_type: str = "evaluation"
    description: str = ""


class ExperimentContext:
    def __init__(self, exp_id: Optional[str] = None) -> None:
        self._data: dict[str, Any] = ExperimentContext.default_values(exp_id)
        self._frozen = False

    @staticmethod
    def default_values(exp_id: Optional[str]) -> Dict[str, Any]:
        from ttex.log import LOGGER_NAME as LOGGER_NAME

        exp_id = exp_id if exp_id is not None else str(uuid4())

        default_dict = {
            "logger_name": LOGGER_NAME,
            "wandb_logger_name": None,
            "coco_logger_name": None,
            "exp_id": exp_id,
            "approach_name": None,
            "artifacts": [],
        }
        return default_dict

    def set(self, key: str, value: Any) -> None:
        if self._frozen:
            raise RuntimeError("Cannot modify frozen context")
        self._data[key] = value

    def get(self, key: str, default: Any = None) -> Any:
        return self._data.get(key, default)

    def add_artifact(self, artifact: Artifact) -> None:
        if "artifacts" not in self._data:
            self._data["artifacts"] = []
        self._data["artifacts"].append(artifact)

    def freeze(self) -> None:
        self._frozen = True
