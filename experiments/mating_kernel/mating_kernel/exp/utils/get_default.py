import inspect


def get_default(cls, param):
    parameter = inspect.signature(cls).parameters.get(param)
    return parameter.default if parameter else None
