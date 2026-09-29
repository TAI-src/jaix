from mating_kernel.pymoo.do_recorder import DoRecorderMixin

import inspect


def make_recorded(operator):
    cls = type(operator)

    if issubclass(cls, DoRecorderMixin):
        return operator

    recorded_cls = type(f"Recorded{cls.__name__}", (DoRecorderMixin, cls), {})

    operator.__class__ = recorded_cls
    return operator


def make_recordable(cls):

    signature = inspect.signature(cls.__init__)

    class RecordableObject(cls):

        def __init__(self, *args, record=None, **kwargs):
            record = set(record or [])

            # Bind arguments exactly as Python would bind them to __init__
            bound = signature.bind(self, *args, **kwargs)
            bound.apply_defaults()

            # Replace requested operators with recording versions
            for name in record:
                if name not in bound.arguments:
                    raise ValueError(
                        f"{name!r} is not an argument of " f"{cls.__name__}.__init__"
                    )

                operator = bound.arguments[name]

                if operator is not None:
                    bound.arguments[name] = make_recorded(operator)

            # Call the original constructor with the modified arguments
            bound.arguments.pop("self")
            super().__init__(**bound.arguments)

    RecordableObject.__name__ = f"Recorded{cls.__name__}"

    return RecordableObject
