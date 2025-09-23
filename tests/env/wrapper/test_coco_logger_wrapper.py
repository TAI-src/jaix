from jaix.env.wrapper import (
    COCOLoggerWrapper,
    COCOLoggerWrapperConfig,
    WrappedEnvFactory as WEF,
)
from . import DummyEnv, TestHandler
import pytest
import shutil
import os.path as osp
from cocopp.pproc import DictAlg

algo_name = "test_algo"


@pytest.mark.parametrize("wef", [True, False])
def test_basic(wef):
    config = COCOLoggerWrapperConfig(algo_name=algo_name)
    config.setup()
    assert config.passthrough
    env = DummyEnv()

    if wef:
        wrapped_env = WEF.wrap(env, [(COCOLoggerWrapper, config)])
    else:
        wrapped_env = COCOLoggerWrapper(config, env)
    assert getattr(wrapped_env, "coco_logger", None) is not None

    test_handler = TestHandler(level="INFO")
    wrapped_env.coco_logger.addHandler(test_handler)

    assert not wrapped_env._started
    wrapped_env.reset()
    msg = test_handler.last_record.getMessage()
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
    assert "data_1/f1_d3_i1.tdat" in msg

    success = config.teardown()
    assert success
    res = config.res
    assert isinstance(res, DictAlg)
    result_dict = res[(algo_name, "")]
    # Check that pproc ran successfully
    assert len(result_dict) > 0

    assert osp.exists(osp.join(config.exp_id, algo_name))
    assert osp.exists(osp.join(config.exp_id, "ppdata"))
    shutil.rmtree(config.exp_id)
