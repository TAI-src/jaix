import inspect

from mating_kernel.pymoo.do_recorder import DoRecorderMixin, RecordingConfig


def make_recorded(operator, recording_config: RecordingConfig):
    cls = type(operator)

    if issubclass(cls, DoRecorderMixin):
        return operator

    recorded_cls = type(f"Recorded{cls.__name__}", (DoRecorderMixin, cls), {})

    operator.__class__ = recorded_cls
    DoRecorderMixin.__init__(operator, recording_config=recording_config)
    return operator


def make_recordable(cls):

    signature = inspect.signature(cls.__init__)

    class RecordableObject(cls):

        def __init__(
            self,
            *args,
            record_args: dict[str, RecordingConfig] | None = None,
            record_attributes: dict[str, RecordingConfig] | None = None,
            **kwargs,
        ):
            record_args = record_args or {}
            record_attributes = record_attributes or {}

            # Bind arguments exactly as Python would bind them to __init__
            bound = signature.bind(self, *args, **kwargs)
            bound.apply_defaults()

            # Replace requested operators with recording versions
            for name, rec_config in record_args.items():
                operator = bound.arguments[name]

                if operator is not None:
                    bound.arguments[name] = make_recorded(
                        operator, recording_config=rec_config
                    )
                    assert isinstance(
                        bound.arguments[name], DoRecorderMixin
                    ), f"{name} is not a DoRecorderMixin"

            # Call the original constructor with the modified arguments
            bound.arguments.pop("self")
            super().__init__(**bound.arguments)

            # Check for any missing arguments that were requested to be recorded
            for name, rec_config in record_attributes.items():
                attr_value = getattr(self, name, None)
                if attr_value is None:
                    raise ValueError(
                        f"Cannot record '{name}': attribute not found in instance."
                    )
                setattr(
                    self,
                    name,
                    make_recorded(attr_value, recording_config=rec_config),
                )
                assert isinstance(
                    getattr(self, name), DoRecorderMixin
                ), f"{name} is not a DoRecorderMixin"

    RecordableObject.__name__ = f"Recorded{cls.__name__}"

    return RecordableObject
