from jaix.env.wrapper.coco_logger_wrapper import (
    COCOLoggerWrapper,
    COCOLoggerWrapperConfig,
)
from jaix.env.wrapper.wrapped_env_factory import WrappedEnvFactory as WEF
from . import DummyEnv, TestHandler
import pytest
import shutil
import os.path as osp
import os
import logging
import jaix.utils.globals as globals
from jaix.utils.exp_id import set_exp_id
from uuid import uuid4

algo_name = "test_algo"


@pytest.mark.parametrize("wef", [True, False])
def test_basic(wef):
    config = COCOLoggerWrapperConfig(algo_name=algo_name)
    config.setup()
    assert config.passthrough

    # Check that logger exists
    assert globals.COCO_LOGGER_NAME in logging.Logger.manager.loggerDict
    # Add test handler to be able to read output
    logger = logging.getLogger(globals.COCO_LOGGER_NAME)
    test_handler = TestHandler(level="INFO")
    logger.addHandler(test_handler)

    # simulate experiment setting experiment id
    exp_id = f"test_coco_logger_wrapper_{uuid4()}"
    set_exp_id(exp_id)

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

    assert osp.exists(osp.join(wrapped_env.exp_id, "DummyEnv", algo_name))
    assert osp.exists(
        osp.join(wrapped_env.exp_id, "DummyEnv", algo_name, "data_1", "f1_d18_i1.tdat")
    )
    assert osp.exists(
        osp.join(wrapped_env.exp_id, "DummyEnv", algo_name, "data_1", "f1_d18_i1.dat")
    )
    assert osp.exists(osp.join(wrapped_env.exp_id, "DummyEnv", algo_name, "f1_i1.info"))
    # assert osp.exists(osp.join(config.exp_id, "ppdata"))

    # clean up folder (needs to remove exdata)
    experiment_dir = wrapped_env.exp_id.split(os.sep)[0]
    shutil.rmtree(experiment_dir)
