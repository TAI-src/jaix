from typing import Optional

_approach_name: Optional[str] = None  # private global variable


def set_approach_name(approach_name: Optional[str]) -> None:
    global _approach_name
    _approach_name = approach_name


def get_approach_name() -> Optional[str]:
    return _approach_name
