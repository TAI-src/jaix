from mating_kernel.exp.utils.get_default import get_default


class DummyClass:
    def __init__(self, param1=10, param2="default", param3=None):
        self.param1 = param1
        self.param2 = param2
        self.param3 = param3


def test_get_default():
    assert get_default(DummyClass, "param1") == 10
    assert get_default(DummyClass, "param2") == "default"
    assert get_default(DummyClass, "param3") is None
    assert get_default(DummyClass, "non_existent_param") is None
