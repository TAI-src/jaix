from typing import Optional

_experiment_id: Optional[str] = None  # private global variable


def set_exp_id(exp_id: Optional[str]) -> None:
    global _experiment_id
    _experiment_id = exp_id


def get_exp_id() -> Optional[str]:
    return _experiment_id
