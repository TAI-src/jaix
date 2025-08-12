from enum import Enum

TaskType = Enum(
    "TaskType", [("C1", "binary"), ("R", "regression"), ("CM", "multiclass")]
)

try:
    from jaix.env.utils.hpo.tabrepo_adapter import TabrepoAdapter
except ImportError:
    # If the import fails, we set TabrepoAdapter to None
    TabrepoAdapter = None  # type: ignore[assignment,misc]
