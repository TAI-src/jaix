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
import logging

algo_name = "test_algo"


@pytest.mark.parametrize("wef", [True, False])
def test_basic(wef):
    config = COCOLoggerWrapperConfig(algo_name=algo_name)
    config.setup()
    assert config.passthrough

    logger = logging.getLogger("coco_logger")
    test_handler = TestHandler(level="INFO")
    logger.addHandler(test_handler)

    env = DummyEnv(dimension=18)
    if wef:
        wrapped_env = WEF.wrap(env, [(COCOLoggerWrapper, config)])
    else:
        wrapped_env = COCOLoggerWrapper(config, env)
    assert getattr(wrapped_env, "coco_logger", None) is not None

    msg = test_handler.last_record.getMessage()
    assert "% f evaluations" in msg
    assert len(test_handler.record_log) == 1

    # Check that reset does not emit new starts
    wrapped_env.reset()
    assert len(test_handler.record_log) == 1

    wrapped_env.step(env.action_space.sample())
    msg = test_handler.last_record.getMessage()
    assert "1 0 " in msg

    wrapped_env.step(env.action_space.sample())
    msg = test_handler.last_record.getMessage()
    assert "2 0 " in msg

    wrapped_env.close()
    msg = test_handler.last_record.getMessage()
    assert "data_1/f1_d18_i1.tdat" in msg

    success = config.teardown()
    assert success
    # TODO: cocopp
    """
    res = config.res
    assert isinstance(res, DictAlg)
    result_dict = res[(algo_name, "")]
    # Check that pproc ran successfully
    assert len(result_dict) > 0
    """

    assert osp.exists(osp.join(config.exp_id, "DummyEnv", algo_name))
    assert osp.exists(
        osp.join(config.exp_id, "DummyEnv", algo_name, "data_1", "f1_d18_i1.tdat")
    )
    assert osp.exists(
        osp.join(config.exp_id, "DummyEnv", algo_name, "data_1", "f1_d18_i1.dat")
    )
    assert osp.exists(osp.join(config.exp_id, "DummyEnv", algo_name, "f1_i1.info"))
    # assert osp.exists(osp.join(config.exp_id, "ppdata"))
    shutil.rmtree(config.exp_id)
