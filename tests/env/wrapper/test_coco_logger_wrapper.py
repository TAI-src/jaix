from jaix.env.wrapper import (
    COCOLoggerWrapper,
    COCOLoggerWrapperConfig,
    WrappedEnvFactory as WEF,
)
from . import DummyEnv, TestHandler
import pytest


@pytest.mark.parametrize("wef", [True, False])
def test_basic(wef):
    config = COCOLoggerWrapperConfig(algo_name="test_algo")
    assert config.passthrough
    env = DummyEnv()

    if wef:
        wrapped_env = WEF.wrap(env, [(COCOLoggerWrapper, config)])
    else:
        wrapped_env = COCOLoggerWrapper(config, env)
    wrapped_env.exp_id = "exp_id"  # TODO: how to set this?
    assert hasattr(wrapped_env, "coco_logger")
    handlers = [
        h for h in wrapped_env.coco_logger.handlers if isinstance(h, TestHandler)
    ]
    if handlers:
        test_handler = handlers[0]
    else:
        test_handler = TestHandler(level="INFO")
        wrapped_env.coco_logger.addHandler(test_handler)
    print(wrapped_env.coco_logger)
    print(wrapped_env.coco_logger.handlers)

    assert not wrapped_env._started
    wrapped_env.reset()
    msg = test_handler.last_record.getMessage()
    print(test_handler.record_log)
    assert "% f evaluations" in msg
    assert wrapped_env._started

    wrapped_env.step(env.action_space.sample())
    msg = test_handler.last_record.getMessage()
    assert "1 0 " in msg

    wrapped_env.step(env.action_space.sample())
    msg = test_handler.last_record.getMessage()
    assert "2 0 " in msg

    wrapped_env.close()
    msg = test_handler.last_record.getMessage()
    assert "data_0/exp_id_0_d3_i1.tdat" in msg
