from jaix.env.wrapper import PassthroughWrapper
import gymnasium as gym
from ttex.config import ConfigurableObject, Config
from typing import List, Optional
from ttex.log.coco import (
    setup_coco_logger,
    COCOStart,
    COCOEnd,
    COCOEval,
    teardown_coco_logger,
)
import logging
import cocopp
from uuid import uuid4
import os.path as osp
import os
import contextlib
from jaix import LOGGER_NAME
import numpy as np

logger = logging.getLogger(LOGGER_NAME)


class COCOLoggerWrapperConfig(Config):
    def __init__(
        self,
        algo_name: str,
        algo_info: str = "",
        exp_id: Optional[str] = None,
        logger_name: str = "coco_logger",
        base_evaluation_triggers: Optional[List[int]] = None,
        number_evaluation_triggers: int = 20,
        improvement_steps: float = 1e-5,
        number_target_triggers: int = 20,
        target_precision: float = 1e-8,
        passthrough: bool = True,
    ):
        self.algo_name = algo_name
        self.exp_id = exp_id if exp_id is not None else str(uuid4())
        self.algo_info = algo_info
        # TODO: potentially add some env info here too
        self.logger_name = logger_name
        self.passthrough = passthrough
        self.base_evaluation_triggers = base_evaluation_triggers
        self.number_evaluation_triggers = number_evaluation_triggers
        self.improvement_steps = improvement_steps
        self.number_target_triggers = number_target_triggers
        self.target_precision = target_precision

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


class COCOLoggerWrapper(PassthroughWrapper, ConfigurableObject):
    config_class = COCOLoggerWrapperConfig

    def __init__(
        self,
        config: COCOLoggerWrapperConfig,
        env: gym.Env,
    ):
        ConfigurableObject.__init__(self, config)
        PassthroughWrapper.__init__(self, env, self.passthrough)
        self.coco_logger = logging.getLogger(self.logger_name)
        self.emit_start()  # Emit start on init

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
            exp_id=self.exp_id,  # Get from wandb
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
        # COCO logger eval
        raw_r = info["raw_r"] if "raw_r" in info else r
        coco_eval = COCOEval(
            x=action,
            mf=raw_r,  # TODO raw_r or observation?
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
