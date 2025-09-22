from jaix.env.wrapper import PassthroughWrapper
import gymnasium as gym
from ttex.config import ConfigurableObject, Config
from typing import List, Optional
from ttex.log.coco import setup_coco_logger, COCOStart, COCOEnd, COCOEval


class COCOLoggerWrapperConfig(Config):
    def __init__(
        self,
        algo_name: str,
        algo_info: str = "",
        logger_name: str = "coco_logger",
        base_evaluation_triggers: Optional[List[int]] = None,
        number_evaluation_triggers: int = 20,
        improvement_steps: float = 1e-5,
        number_target_triggers: int = 20,
        target_precision: float = 1e-8,
        passthrough: bool = True,
    ):
        self.algo_name = algo_name
        self.algo_info = algo_info
        self.logger_name = logger_name
        self.passthrough = passthrough
        self.base_evaluation_triggers = base_evaluation_triggers
        self.number_evaluation_triggers = number_evaluation_triggers
        self.improvement_steps = improvement_steps
        self.number_target_triggers = number_target_triggers
        self.target_precision = target_precision


class COCOLoggerWrapper(PassthroughWrapper, ConfigurableObject):
    config_class = COCOLoggerWrapperConfig

    def __init__(
        self,
        config: COCOLoggerWrapperConfig,
        env: gym.Env,
    ):
        ConfigurableObject.__init__(self, config)
        PassthroughWrapper.__init__(self, env, self.passthrough)
        self.coco_logger = setup_coco_logger(
            name=self.logger_name,
            base_evaluation_triggers=self.base_evaluation_triggers,
            number_evaluation_triggers=self.number_evaluation_triggers,
            improvement_steps=self.improvement_steps,
            number_target_triggers=self.number_target_triggers,
            target_precision=self.target_precision,
        )
        self._exp_id = None
        self._started = False

    @property
    def exp_id(self):
        return self._exp_id

    @exp_id.setter
    def exp_id(self, value: str):
        self._exp_id = value

    def reset(self, **kwargs):
        obs, info = self.env.reset(**kwargs)
        if self._started:  # If previously started, need to restart
            # Add COCO end
            self.coco_logger.info(COCOEnd())
        self.emit_start()
        return obs, info

    def emit_start(self):
        # Tell COCO that a new experiment is starting
        coco_start = COCOStart(
            algo=self.algo_name,
            problem=(
                self.env.unwrapped.func_id
                if hasattr(self.env.unwrapped, "func_id")
                else 0
            ),
            dim=self.action_space.shape[0],
            inst=self.env.unwrapped.inst if hasattr(self.env.unwrapped, "inst") else 0,
            suite=self.env.unwrapped.__class__.__name__,
            exp_id=self.exp_id,
            algo_info=self.algo_info,
            fopt=(
                self.env.unwrapped.fopt if hasattr(self.env.unwrapped, "fopt") else None
            ),
        )
        print(
            f"algo: {self.algo_name}, problem: {coco_start.problem}, dim: {coco_start.dim}, inst: {coco_start.inst}, suite: {coco_start.suite}"
        )

        self.coco_logger.info(coco_start)
        self._started = True
        return coco_start

    def step(self, action):
        (
            obs,
            r,
            term,
            trunc,
            info,
        ) = self.env.step(action)
        if not self._started:
            self.emit_start()
        # COCO logger eval
        raw_r = info["raw_r"] if "raw_r" in info else r
        coco_eval = COCOEval(
            x=action,
            mf=raw_r,  # TODO raw_r or observation?
        )
        self.coco_logger.info(coco_eval)
        return obs, r, term, trunc, info

    def close(self):
        self.env.close()
        # Tell COCO that the experiment is done
        self.coco_logger.info(COCOEnd())
        self._started = False
