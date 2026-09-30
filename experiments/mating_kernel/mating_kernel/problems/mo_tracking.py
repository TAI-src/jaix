from typing import Any, cast

import numpy as np
from jaix.env.singular.ec_env import ECEnvironment, ECEnvironmentConfig
from jaix.env.utils.archive.entry_scorer import (
    ReferenceVectorDistanceScorer,
)
from jaix.env.utils.archive.mo_archive import (
    KeepDominated,
    MOArchive,
    MOArchiveConfig,
    MOArchiveEntry,
)
from jaix.env.utils.problem.static_problem import StaticProblem
from typing import TypeVar


class MOEvalEntry(MOArchiveEntry):
    def __init__(self, x: np.ndarray, y: np.ndarray):
        self.x = x
        self.y = y

    def parse(self) -> np.ndarray:
        return self.y


class MOTrackingMixin:

    def __init_subclass__(cls, **kwargs):
        """
        This makes sure that the next class in the MRO is a subclass of StaticProblem.
        """
        super().__init_subclass__(**kwargs)

        mro = cls.mro()
        mixin_index = mro.index(MOTrackingMixin)

        if mixin_index + 1 >= len(mro):
            raise TypeError("MOTrackingMixin must be followed by another class")

        next_class = mro[mixin_index + 1]

        if not issubclass(next_class, StaticProblem):
            raise TypeError(
                f"{cls.__name__}: MOTrackingMixin must be followed by "
                f"ExpectedClass, got {next_class.__name__}"
            )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.archive = self._create_eval_archive(self)

    @staticmethod
    def _create_eval_archive(func: StaticProblem) -> MOArchive:
        archive_config = MOArchiveConfig(
            MOEvalEntry,
            secondary_criterion_class=ReferenceVectorDistanceScorer,
            max_size=None,
            keep_dominated=KeepDominated.NONE,
            only_new_entries=False,
            num_refpoints="original",
        )
        env = ECEnvironment(ECEnvironmentConfig(budget_multiplier=1), func=func)
        return MOArchive(archive_config, env=env)

    def _eval(self, x) -> tuple[list[float], list[float]]:
        """
        Evaluate the objective function.
            :param x: The input vector.
            :return: Tuple of objective function value (clean and noisy).
        """
        y_raw, y_noisy = cast(Any, super())._eval(x)
        entry = MOEvalEntry(x=np.array(x), y=np.array(y_raw))
        self.archive.add([entry], queue=True)
        return y_raw, y_noisy

    def get_archive_stats(self) -> dict[str, Any]:
        self.archive.add_queued()
        stats = self.archive.get_archive_stats()
        stats["fevals"] = getattr(self, "evaluations", None)
        return stats


T = TypeVar("T", bound=StaticProblem)


def make_tracked(cls: type[T]) -> type[T]:
    if not issubclass(cls, StaticProblem):
        raise TypeError(f"make_tracked expects a StaticProblem subclass, got {cls!r}")

    return type(
        f"Tracked{cls.__name__}",
        (MOTrackingMixin, cls),
        {
            "__module__": cls.__module__,
            "__doc__": cls.__doc__,
        },
    )
