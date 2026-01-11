from jaix.utils.experiment_context import ExperimentContext, Artifact


def test_artifact():
    art = Artifact(
        name="test_artifact",
        local_path="/tmp/test_artifact",
        artifact_type="dataset",
        description="This is a test artifact",
    )
    assert art.name == "test_artifact"
    assert art.local_path == "/tmp/test_artifact"
    assert art.artifact_type == "dataset"
    assert art.description == "This is a test artifact"


def test_default_values():
    def_dict = ExperimentContext.default_values(None)
    assert "exp_id" in def_dict
    assert "logger_name" in def_dict
    assert def_dict["exp_id"] is not None
    assert def_dict["logger_name"] is not None
    assert def_dict["wandb_logger_name"] is None
    assert def_dict["coco_logger_name"] is None
    assert def_dict["approach_name"] is None
    assert def_dict["artifacts"] == []

    def_dict = ExperimentContext.default_values("custom_exp_id")
    assert def_dict["exp_id"] == "custom_exp_id"


def test_set_get():
    ctx = ExperimentContext()
    ctx.set("test_key", "test_value")
    assert ctx.get("test_key") == "test_value"
    assert ctx.get("non_existent_key", "default_value") == "default_value"


def test_freeze():
    ctx = ExperimentContext()
    ctx.set("key1", "value1")
    ctx.freeze()
    try:
        ctx.set("key2", "value2")
        assert False, "Expected RuntimeError when modifying frozen context"
    except RuntimeError:
        pass  # Expected


def test_add_artifact():
    ctx = ExperimentContext()
    art1 = Artifact(name="artifact1", local_path="/path/to/artifact1")
    art2 = Artifact(name="artifact2", local_path="/path/to/artifact2")
    ctx.add_artifact(art1)
    ctx.freeze()  # Freezing should not prevent adding artifacts
    ctx.add_artifact(art2)
    artifacts = ctx.get("artifacts")
    assert len(artifacts) == 2
    assert artifacts[0].name == "artifact1"
    assert artifacts[1].name == "artifact2"
