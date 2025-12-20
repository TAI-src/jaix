from jaix.env.wrapper.passthrough_wrapper import PassthroughWrapper
import gymnasium as gym
from ttex.config import ConfigurableObject, Config
from typing import List, Optional
from ttex.log.coco import (
    COCOStart,
    COCOEnd,
    COCOEval,
)
from ttex.log import setup_coco_logger, teardown_coco_logger
import logging
import os.path as osp
import numpy as np
from jaix.utils.exp_id import get_exp_id
from jaix.utils.approach_name import get_approach_name
from jaix.env.wrapper.value_track_wrapper import ValueTrackWrapper

import jaix.utils.globals as globals

logger = logging.getLogger(globals.LOGGER_NAME)


class COCOLoggerWrapperConfig(Config):
    def __init__(
        self,
        algo_name: Optional[str] = None,
        algo_info: str = "",
        logger_name: Optional[str] = None,
        base_evaluation_triggers: Optional[List[int]] = None,
        number_evaluation_triggers: int = 20,
        improvement_steps: float = 1e-5,
        number_target_triggers: int = 20,
        target_precision: float = 1e-8,
        passthrough: bool = True,
        state_eval: str = "obs0",  # Which value should be logged
        is_min: bool = True,  # Whether lower is better for state_eval
    ):
        self.algo_name = algo_name if algo_name is not None else get_approach_name()
        self.algo_info = algo_info
        # TODO: potentially add some env info here too
        self.logger_name = (
            logger_name if logger_name is not None else globals.COCO_LOGGER_NAME
        )
        if self.logger_name == globals.LOGGER_NAME:
            raise ValueError(
                "COCOLoggerWrapperConfig: logger_name cannot be the root logger name."
            )
        globals.COCO_LOGGER_NAME = self.logger_name
        self.passthrough = passthrough
        self.base_evaluation_triggers = base_evaluation_triggers
        self.number_evaluation_triggers = number_evaluation_triggers
        self.improvement_steps = improvement_steps
        self.number_target_triggers = number_target_triggers
        self.target_precision = target_precision
        self.state_eval = state_eval
        self.is_min = is_min

    def _setup(self):
        setup_coco_logger(
            name=self.logger_name,
            base_evaluation_triggers=self.base_evaluation_triggers,
            number_evaluation_triggers=self.number_evaluation_triggers,
            improvement_steps=self.improvement_steps,
            number_target_triggers=self.number_target_triggers,
            target_precision=self.target_precision,
        )
        return True

    def _teardown(self):
        # This also triggers writing the files
        teardown_coco_logger(self.logger_name)

        # If results are generated, run cocopp post-processing
        # TODO: set up cocopp
        """
        if osp.exists(osp.join(self.exp_id, self.algo_name)):
            # Run cocopp post-processing on the generated files (but quietly)
            with open(os.devnull, "w") as devnull:
                with contextlib.redirect_stdout(devnull):
                    self.res = cocopp.main(
                        f"-o {osp.join(self.exp_id, 'ppdata')} {osp.join(self.exp_id, self.algo_name)}"
                    )
        else:
            logger.warning(
                f"No results found in {osp.join(self.exp_id, self.algo_name)}. Skipping cocopp post-processing."
            )
            self.res = None
        """
        return True


class COCOLoggerWrapper(ConfigurableObject, ValueTrackWrapper):
    """
    A wrapper that logs environment interactions using COCO format for later postprocessing.
    """

    config_class = COCOLoggerWrapperConfig

    def __init__(
        self,
        config: COCOLoggerWrapperConfig,
        env: gym.Env,
    ):
        ConfigurableObject.__init__(self, config)
        ValueTrackWrapper.__init__(
            self,
            env,
            passthrough=config.passthrough,
            state_eval=config.state_eval,
            is_min=config.is_min,
        )
        self.coco_logger = logging.getLogger(self.logger_name)
        exp_id = get_exp_id()
        assert exp_id is not None, "COCOLoggerWrapper: exp_id must be set globally"
        self.exp_id = COCOLoggerWrapper.coco_dir(exp_id)
        self.emit_start()  # Emit start on init

    @staticmethod
    def coco_dir(exp_id: Optional[str]) -> str:
        exp_id = exp_id if exp_id is not None else get_exp_id()
        assert exp_id is not None, "exp_id must be provided or set globally"
        return osp.join(exp_id, "coco_exdata")

    def emit_start(self):
        # Tell COCO that a new experiment is starting
        constant_dim = not hasattr(self.env, "constant_dim") or self.env.constant_dim
        suite_name = (  # Especially important for composite envs
            self.suite_name
            if hasattr(self, "suite_name")
            else type(self.env.unwrapped).__name__
        )
        coco_start = COCOStart(
            algo=self.algo_name,
            problem=(
                self.env.unwrapped.func_id + 1  # TODO: fix for 0-indexing
                if hasattr(self.env.unwrapped, "func_id")
                else 1
            ),
            dim=np.prod(self.action_space.shape) if constant_dim else 0,
            inst=self.env.unwrapped.inst if hasattr(self.env.unwrapped, "inst") else 1,
            suite=suite_name,
            exp_id=self.exp_id,
            algo_info=self.algo_info,
            fopt=(
                self.env.unwrapped.fopt if hasattr(self.env.unwrapped, "fopt") else None
            ),
        )
        self.coco_logger.info(coco_start)
        logger.debug(f"COCOStart emitted: {coco_start} {self.exp_id}")
        return coco_start

    def step(self, action):
        (
            obs,
            r,
            term,
            trunc,
            info,
        ) = self.env.step(action)

        # Get the value for tracking and update internal state
        val = self.get_val(obs, r, info, self.state_eval)
        self.update_vals(val)

        coco_eval = COCOEval(
            x=action,
            mf=val,  # TODO: should also be logging noisy values
        )
        self.coco_logger.info(coco_eval)
        logger.debug(f"COCOEval emitted: {coco_eval} {self.exp_id}")
        return obs, r, term, trunc, info

    def close(self):
        self.env.close()
        # Tell COCO that the experiment is done
        self.coco_logger.info(COCOEnd())
        logger.debug(f"COCOEnd emitted {self.exp_id}")
        return True
