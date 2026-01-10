from jaix.env.wrapper.wandb_wrapper import WandbWrapper, WandbWrapperConfig
from jaix.env.wrapper.wrapped_env_factory import (
    WrappedEnvFactory as WEF,
)
from jaix.env.wrapper.improvement_reward_wrapper import (
    ImprovementRewardWrapper,
    ImprovementRewardWrapperConfig,
)
from . import DummyEnv, DummyWrapper, DummyWrapperConfig, TestHandler
from gymnasium.utils.env_checker import check_env
import ast
import pytest
import logging
from ttex.log import teardown_wandb_logger
from jaix.utils.experiment_context import ExperimentContext
from jaix.utils.wandb_session import WandbSession
import os
from shutil import rmtree

WANDB_LOGGER_TEST_NAME = "WandbLoggerTest"


@pytest.fixture(autouse=True)
def run_around_tests():
    # ensure wandb offline mode
    prev_mode = os.environ.get("WANDB_MODE", "online")
    os.environ["WANDB_MODE"] = "offline"
    yield
    # Code that will run after your test, e.g. teardown
    os.environ["WANDB_MODE"] = prev_mode


@pytest.fixture(scope="module", autouse=True)
def after_all_tests():
    # remove wandb logs
    # setup (optional)
    yield
    # this runs once after all tests in the file
    teardown_wandb_logger(WANDB_LOGGER_TEST_NAME)
    rmtree("wandb", ignore_errors=True)


def prep_config(**kwargs):
    config = WandbWrapperConfig(wandb_logger_name=WANDB_LOGGER_TEST_NAME, **kwargs)
    ctx = ExperimentContext()
    config.setup(ctx)
    wandb_session = WandbSession(config)
    wandb_session.start(ctx)
    ctx.freeze()
    config.set_context(ctx)
    test_handler = TestHandler(level="INFO")
    logger = logging.getLogger(WANDB_LOGGER_TEST_NAME)
    logger.addHandler(test_handler)
    return config, test_handler


@pytest.mark.parametrize("wef", [True, False])
def test_basic(wef):
    config, test_handler = prep_config()
    env = DummyEnv()

    if wef:
        wrapped_env = WEF.wrap(env, [(WandbWrapper, config)])
    else:
        wrapped_env = WandbWrapper(config, env)
    assert hasattr(wrapped_env, "wandb_logger")

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
    double_name = "test_logger"
    ctx = ExperimentContext()
    ctx.set("logger_name", double_name)
    config = WandbWrapperConfig(wandb_logger_name=double_name)
    with pytest.raises(ValueError):
        config.setup(ctx)


def test_additions():
    config, test_handler = prep_config()
    env = ImprovementRewardWrapper(
        ImprovementRewardWrapperConfig(state_eval="obs0"), DummyEnv()
    )  # Adds raw_r
    env = DummyWrapper(DummyWrapperConfig(), env)  # Adds env_step
    wrapped_env = WandbWrapper(config, env)

    wrapped_env.reset()
    wrapped_env.step(wrapped_env.action_space.sample())

    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert msg["env/step"] == msg["env/log_step"] + 1
    assert "env/raw_obs0/DummyEnv/0/1" in msg
    assert "env/best_raw_obs0/DummyEnv/0/1" in msg


def test_close():
    config, test_handler = prep_config()
    env = DummyEnv()
    wrapped_env = WandbWrapper(config, env)

    wrapped_env.reset()
    wrapped_env.step(wrapped_env.action_space.sample())

    wrapped_env.close()
    msg = ast.literal_eval(test_handler.last_record.getMessage())
    assert "env/close/DummyEnv/0/1/funcs" in msg


# Hacky, but needs to happen before the config test with teardown
@pytest.mark.parametrize(
    "state_eval_imp, state_eval_wandb",
    [("obs0", "obs0"), ("obs0", "r"), ("val", "obs0"), ("val", "r")],
)
def test_wandb_improvement_interaction_v1(state_eval_imp, state_eval_wandb):
    env = ImprovementRewardWrapper(
        ImprovementRewardWrapperConfig(state_eval=state_eval_imp, is_min=True),
        DummyEnv(),
    )
    config, test_handler = prep_config(state_eval=state_eval_wandb, is_min=True)
    wrapped_env = WandbWrapper(config, env)  # Adds raw_obs0

    assert hasattr(wrapped_env, "wandb_logger")

    wrapped_env.reset()
    for _ in range(3):
        wrapped_env.step(wrapped_env.action_space.sample())
    obs, r, term, trunc, info = wrapped_env.step(wrapped_env.action_space.sample())

    gmsg = test_handler.last_record.getMessage()
    msg = ast.literal_eval(gmsg)
    print(msg)

    if state_eval_imp == state_eval_wandb:
        assert wrapped_env.last_val == wrapped_env.env.last_val
    assert info[f"best_raw_{state_eval_imp}"] == wrapped_env.env.best_val
    assert info[f"raw_{state_eval_imp}"] == wrapped_env.env.last_val

    assert msg[f"env/best_raw_{state_eval_wandb}/DummyEnv/0/1"] == wrapped_env.best_val
    assert msg[f"env/raw_{state_eval_wandb}/DummyEnv/0/1"] == wrapped_env.last_val

    assert "env/r/DummyEnv/0/1" in msg

    assert r is not None  # should be the improvement reward
    assert "improvement" in info


@pytest.mark.parametrize(
    "state_eval_imp, state_eval_wandb",
    [("obs0", "obs0"), ("obs0", "r"), ("r", "obs0"), ("r", "r")],
)
def test_wandb_improvement_interaction_v2(state_eval_imp, state_eval_wandb):
    env = DummyEnv()
    config, test_handler = prep_config(state_eval=state_eval_wandb, is_min=True)
    wrapped_env = WandbWrapper(config, env)  # Adds raw_obs0
    assert hasattr(wrapped_env, "wandb_logger")

    # Now wrap with ImprovementRewardWrapper
    imp_config = ImprovementRewardWrapperConfig(state_eval=state_eval_imp, is_min=True)
    wrapped_env = ImprovementRewardWrapper(imp_config, wrapped_env)

    wrapped_env.reset()
    for _ in range(3):
        wrapped_env.step(wrapped_env.action_space.sample())
    obs, r, term, trunc, info = wrapped_env.step(wrapped_env.action_space.sample())

    gmsg = test_handler.last_record.getMessage()
    msg = ast.literal_eval(gmsg)
    assert "env/r/DummyEnv/0/1" in msg
    assert msg["env/r/DummyEnv/0/1"] != r

    assert r is not None  # should be the improvement reward
    assert "improvement" in info

    if state_eval_imp == state_eval_wandb:
        assert wrapped_env.last_val == wrapped_env.env.last_val
    assert info[f"best_raw_{state_eval_imp}"] == wrapped_env.best_val
    assert info[f"raw_{state_eval_imp}"] == wrapped_env.last_val

    assert (
        msg[f"env/best_raw_{state_eval_wandb}/DummyEnv/0/1"] == wrapped_env.env.best_val
    )
    assert msg[f"env/raw_{state_eval_wandb}/DummyEnv/0/1"] == wrapped_env.env.last_val


# changes the logger, needs to happen last
def test_wandb_config():
    config = WandbWrapperConfig(
        wandb_logger_name="WandbLogger",
        custom_metrics={"test_metric": 42},
        snapshot=False,
        snapshot_sensitive_keys=["secret"],
        project="test_project",
        group="test_group",
    )

    ctx = ExperimentContext()
    config.setup(ctx)
    env = DummyEnv()
    wrapped_env = WandbWrapper(config, env)
    assert wrapped_env.wandb_logger.name == "WandbLogger"
    assert ctx.get("wandb_logger_name") == "WandbLogger"

    logger = logging.getLogger("WandbLogger")
    assert logger._wandb_setup
    teardown_wandb_logger("WandbLogger")
    assert not logger._wandb_setup
