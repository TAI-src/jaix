from collections.abc import Callable, Iterable
from functools import wraps
from typing import Any, TypeVar

T = TypeVar("T")


def marks_dirty(method: Callable[..., Any]) -> Callable[..., Any]:
    @wraps(method)
    def wrapper(self: "DirtyList[T]", *args: Any, **kwargs: Any) -> Any:
        result = method(self, *args, **kwargs)
        self._dirty = True
        return result

    return wrapper


class DirtyList(list[T]):
    def __init__(self, iterable: Iterable[T] = ()):
        super().__init__(iterable)
        self._dirty = bool(self)  # List is dirty if it has initial elements

    @property
    def dirty(self) -> bool:
        return self._dirty

    def mark_clean(self) -> None:
        self._dirty = False

    @marks_dirty
    def append(self, item: T) -> None:
        super().append(item)

    @marks_dirty
    def extend(self, iterable: Iterable[T]) -> None:
        super().extend(iterable)

    @marks_dirty
    def insert(self, index: int, item: T) -> None:
        super().insert(index, item)

    @marks_dirty
    def remove(self, item: T) -> None:
        super().remove(item)

    @marks_dirty
    def clear(self) -> None:
        super().clear()

    @marks_dirty
    def pop(self, index: int = -1) -> T:
        return super().pop(index)
