from jaix.utils.wandb_session import WandbSession
import pytest
import os
from shutil import rmtree
from jaix.utils.experiment_context import ExperimentContext, Artifact
from ttex.config import Config
from ttex.log import setup_wandb_logger

WANDB_LOGGER_TEST_NAME = "wandb_test_logger"


@pytest.fixture(autouse=True)
def run_around_tests():
    # ensure wandb offline mode
    prev_mode = os.environ.get("WANDB_MODE", "online")
    os.environ["WANDB_MODE"] = "offline"

    setup_wandb_logger(snapshot=False, name=WANDB_LOGGER_TEST_NAME)
    yield
    # Code that will run after your test, e.g. teardown
    os.environ["WANDB_MODE"] = prev_mode


@pytest.fixture(scope="module", autouse=True)
def after_all_tests():
    # remove wandb logs
    # setup (optional)
    yield
    # this runs once after all tests in the file
    rmtree("wandb", ignore_errors=True)


class DummyConfig:
    def to_dict(self):
        return {"param1": 10, "param2": "value"}


def prep():
    config = DummyConfig()
    ctx = ExperimentContext()
    ctx.set("wandb_logger_name", WANDB_LOGGER_TEST_NAME)
    return config, ctx


def test_wandb_session_start_end():
    config, ctx = prep()
    wandb_session = WandbSession(config)
    wandb_session.start(ctx)

    exp_id = ctx.get("exp_id")
    assert exp_id is not None
    assert wandb_session.run.id == exp_id

    assert wandb_session.config_dict == DummyConfig().to_dict()

    # Add dummy artifacts to context
    artifact1 = Artifact(name="artifact1", local_path="/path/to/artifact1")
    artifact2 = Artifact(
        name="artifact2",
        local_path="/path/to/artifact2",
        artifact_type="model",
        description="Test model artifact",
    )
    ctx.add_artifact(artifact1)
    ctx.add_artifact(artifact2)

    wandb_session.end(ctx)
    # Since we are in offline mode, we cannot check wandb server,
    # but we can ensure no exceptions were raised during the process.
    # If needed, more detailed checks can be added with a mock wandb logger.
    # Check that artifacts are still in context after ending session
    artifacts = ctx.get("artifacts", [])
    assert len(artifacts) == 2
    assert artifacts[0].name == "artifact1"
    assert artifacts[1].name == "artifact2"


def test_wandb_session_no_artifacts():
    config, ctx = prep()
    wandb_session = WandbSession(config)
    wandb_session.start(ctx)

    exp_id = ctx.get("exp_id")
    assert exp_id is not None
    assert wandb_session.run.id == exp_id

    # No artifacts added to context
    wandb_session.end(ctx)
    # Ensure no exceptions were raised during the process.
    # If needed, more detailed checks can be added with a mock wandb logger.
    # Check that artifacts list is empty
    artifacts = ctx.get("artifacts", [])
    assert len(artifacts) == 0


def test_robustness():
    # Test with a config that does not implement to_dict
    class EmptyConfig(Config):
        def __init__(self):
            pass

    config = EmptyConfig()
    # Test with experiment context without wandblogger
    ctx = ExperimentContext()
    wandb_session = WandbSession(config)
    wandb_session.start(ctx)

    exp_id = ctx.get("exp_id")
    assert exp_id is not None
    assert wandb_session.run is None

    assert wandb_session.config_dict == {}

    wandb_session.end(ctx)
    # Ensure no exceptions were raised during the process.
    artifacts = ctx.get("artifacts", [])
    assert len(artifacts) == 0
