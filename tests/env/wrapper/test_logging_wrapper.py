from jaix.env.wrapper.wandb_wrapper import WandbWrapper, WandbWrapperConfig
from jaix.env.wrapper.wrapped_env_factory import (
    WrappedEnvFactory as WEF,
)
from jaix.env.wrapper.any_fit_wrapper import AnyFitWrapper
from . import DummyEnv, test_handler, DummyWrapper
from gymnasium.utils.env_checker import check_env
import ast
import pytest


@pytest.mark.parametrize("wef", [True, False])
def test_basic(wef):
    config = WandbWrapperConfig(logger_name="DefaultLogger")
    assert config.passthrough
    env = DummyEnv()

    if wef:
        wrapped_env = WEF.wrap(env, [(WandbWrapper, config)])
    else:
        wrapped_env = WandbWrapper(config, env)
    assert hasattr(wrapped_env, "logger")

    check_env(wrapped_env, skip_render_check=True)

    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert "env/r/DummyEnv/0/1" in msg
    steps = msg["env/step"]
    resets = msg["env/resets/DummyEnv/0/1"]

    wrapped_env.step(wrapped_env.action_space.sample())
    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert msg["env/step"] == steps + 1

    wrapped_env.reset()
    wrapped_env.step(wrapped_env.action_space.sample())
    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert msg["env/resets/DummyEnv/0/1"] == resets + 1


def test_additions():
    config = WandbWrapperConfig(logger_name="DefaultLogger")
    env = AnyFitWrapper(DummyEnv())  # Adds raw_r
    env = DummyWrapper(env)  # Adds env_step
    wrapped_env = WandbWrapper(config, env)

    wrapped_env.reset()
    wrapped_env.step(wrapped_env.action_space.sample())

    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert msg["env/step"] == msg["env/log_step"] + 1
    assert "env/raw_r/DummyEnv/0/1" in msg
    assert "env/best_raw_r/DummyEnv/0/1" in msg


def test_close():
    config = WandbWrapperConfig(logger_name="DefaultLogger")
    env = DummyEnv()
    wrapped_env = WandbWrapper(config, env)

    wrapped_env.reset()
    wrapped_env.step(wrapped_env.action_space.sample())

    wrapped_env.close()
    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert "env/close/DummyEnv/0/1/funcs" in msg
