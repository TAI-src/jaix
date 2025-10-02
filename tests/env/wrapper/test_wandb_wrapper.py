from jaix.env.wrapper.wandb_wrapper import WandbWrapper, WandbWrapperConfig
from jaix.env.wrapper.wrapped_env_factory import (
    WrappedEnvFactory as WEF,
)
from jaix.env.wrapper.any_fit_wrapper import AnyFitWrapper
from . import DummyEnv, test_handler, DummyWrapper, DummyWrapperConfig
from gymnasium.utils.env_checker import check_env
import ast
import pytest
import jaix.utils.globals as globals
import logging
from ttex.log import teardown_wandb_logger


@pytest.fixture(autouse=True)
def run_around_tests():
    prev_logger_name = globals.WANDB_LOGGER_NAME
    globals.WANDB_LOGGER_NAME = globals.LOGGER_NAME
    globals.LOGGER_NAME = "root"  # ensure
    # we use the root logger just for these tests
    yield
    # Code that will run after your test, e.g. teardown
    globals.LOGGER_NAME = globals.WANDB_LOGGER_NAME
    globals.WANDB_LOGGER_NAME = prev_logger_name


@pytest.mark.parametrize("wef", [True, False])
def test_basic(wef):
    config = WandbWrapperConfig()
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


def test_name_conflict():
    with pytest.raises(ValueError):
        config = WandbWrapperConfig(logger_name=globals.LOGGER_NAME)


def test_additions():
    config = WandbWrapperConfig()
    env = AnyFitWrapper(DummyEnv())  # Adds raw_r
    env = DummyWrapper(DummyWrapperConfig(), env)  # Adds env_step
    wrapped_env = WandbWrapper(config, env)

    wrapped_env.reset()
    wrapped_env.step(wrapped_env.action_space.sample())

    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert msg["env/step"] == msg["env/log_step"] + 1
    assert "env/raw_r/DummyEnv/0/1" in msg
    assert "env/best_raw_r/DummyEnv/0/1" in msg


def test_close():
    config = WandbWrapperConfig()
    env = DummyEnv()
    wrapped_env = WandbWrapper(config, env)

    wrapped_env.reset()
    wrapped_env.step(wrapped_env.action_space.sample())

    wrapped_env.close()
    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert "env/close/DummyEnv/0/1/funcs" in msg


def test_wandb_config():
    config = WandbWrapperConfig(
        logger_name="WandbLogger",
        custom_metrics={"test_metric": 42},
        snapshot=False,
        snapshot_sensitive_keys=["secret"],
        project="test_project",
        group="test_group",
    )

    env = DummyEnv()
    wrapped_env = WandbWrapper(config, env)
    assert wrapped_env.logger.name == "WandbLogger"
    assert globals.WANDB_LOGGER_NAME == "WandbLogger"

    config.setup()
    logger = logging.getLogger("WandbLogger")
    assert logger._wandb_setup
    teardown_wandb_logger("WandbLogger")
    assert not logger._wandb_setup
