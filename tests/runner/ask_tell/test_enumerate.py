from jaix.runner.ask_tell.strategy.enumerate import (
    EnumerateATStrat,
    EnumerateATStratConfig,
)
from . import DummyEnv, loop
import pytest
from gymnasium import spaces


def test_init_fail():
    env = DummyEnv()
    # Default dummy env is not MultiBinary
    with pytest.raises(AssertionError):
        EnumerateATStrat(EnumerateATStratConfig(ask_size=5), env)


def test_init():
    env = DummyEnv(
        action_space=spaces.MultiBinary(3),
    )
    enumerate_opt = EnumerateATStrat(EnumerateATStratConfig(ask_size=5), env)
    assert not hasattr(enumerate_opt, "xcurrent")
    assert enumerate_opt.current == 0
    assert enumerate_opt.options == [
        (0, 0, 0),
        (0, 0, 1),
        (0, 1, 0),
        (0, 1, 1),
        (1, 0, 0),
        (1, 0, 1),
        (1, 1, 0),
        (1, 1, 1),
    ]
    assert enumerate_opt.ask_size == 5
    assert env.action_space.contains(enumerate_opt.options[0])


def test_ask():
    env = DummyEnv(
        action_space=spaces.MultiBinary(3),
    )
    enumerate_opt = EnumerateATStrat(EnumerateATStratConfig(ask_size=5), env)
    xcurrent = enumerate_opt.ask(env)
    assert len(xcurrent) == 5  # ask size
    assert len(xcurrent[0]) == 3  # action space.n
    assert enumerate_opt.current == 5  # 5 options asked


def test_ask_end():
    env = DummyEnv(
        action_space=spaces.MultiBinary(3),
    )
    enumerate_opt = EnumerateATStrat(EnumerateATStratConfig(ask_size=5), env)
    # Move to the end
    enumerate_opt.current = len(enumerate_opt.options)
    with pytest.raises(AssertionError):
        enumerate_opt.ask(env)


def test_ask_overflow():
    env = DummyEnv(
        action_space=spaces.MultiBinary(3),
    )
    enumerate_opt = EnumerateATStrat(EnumerateATStratConfig(ask_size=5), env)
    # Move to just before the end
    enumerate_opt.current = len(enumerate_opt.options) - 3
    assert not enumerate_opt.stop()

    xcurrent = enumerate_opt.ask(env)
    assert len(xcurrent) == 3
    assert enumerate_opt.current >= len(enumerate_opt.options)
    assert enumerate_opt.stop()


def test_loop():
    n = 4
    env = DummyEnv(
        action_space=spaces.MultiBinary(4),
    )
    enumerate_opt = EnumerateATStrat(EnumerateATStratConfig(ask_size=5), env)
    sum = 0
    while not enumerate_opt.stop():
        X, _ = loop(None, None, enumerate_opt, env=env)
        sum += len(X)
        assert len(X) <= 5  # ask size
    assert sum == len(enumerate_opt.options)
    assert sum == 2**n  # all combinations of 4 bits
