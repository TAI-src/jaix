from jaix.env.wrapper import LoggingWrapperConfig, LoggingWrapper
from . import DummyEnv, test_handler
from gymnasium.utils.env_checker import check_env
import ast


def test_basic():
    config = LoggingWrapperConfig(logger_name="DefaultLogger")
    env = DummyEnv()
    wrapped_env = LoggingWrapper(config, env)
    assert hasattr(wrapped_env, "logger")

    check_env(wrapped_env, skip_render_check=True)

    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert "env" in msg
    steps = msg["step"]
    assert "reward" in msg

    wrapped_env.step(wrapped_env.action_space.sample())
    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert msg["step"] == steps + 1

    wrapped_env.reset()
    wrapped_env.step(wrapped_env.action_space.sample())
    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert msg["step"] == 1
