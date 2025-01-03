from jaix.env.wrapper import (
    LoggingWrapperConfig,
    LoggingWrapper,
    WrappedEnvFactory as WEF,
)
from . import DummyEnv, test_handler
from gymnasium.utils.env_checker import check_env
import ast
import pytest


@pytest.mark.parametrize("wef", [True, False])
def test_basic(wef):
    config = LoggingWrapperConfig(logger_name="DefaultLogger")
    assert config.passthrough
    env = DummyEnv()

    if wef:
        wrapped_env = WEF.wrap(env, [(LoggingWrapper, config)])
    else:
        wrapped_env = LoggingWrapper(config, env)
    assert hasattr(wrapped_env, "logger")

    check_env(wrapped_env, skip_render_check=True)

    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert "env/r/DummyEnv" in msg
    steps = msg["env/step"]
    resets = msg["env/resets/DummyEnv"]

    wrapped_env.step(wrapped_env.action_space.sample())
    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert msg["env/step"] == steps + 1

    wrapped_env.reset()
    wrapped_env.step(wrapped_env.action_space.sample())
    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert msg["env/resets/DummyEnv"] == resets + 1
