from jaix.env.wrapper import PassthroughWrapper
from . import DummyWrapper, DummyEnv


def test_stopping():
    env = DummyEnv()
    env = DummyWrapper(env)
    assert not env.stop()
    env.stop_dict = {"wrapper 1": "test"}
    env = DummyWrapper(env)
    env = DummyWrapper(env)
    env.stop_dict = {"wrapper 3": "test"}
    assert env.stop()
    assert env.stop() == {"wrapper 1": "test", "wrapper 3": "test"}
    assert not env.unwrapped._stop
    env.unwrapped._stop = True
    assert env.stop()["unwrapped"]
