"""An Experiment"""
from ttex.config import (
    ConfigurableObjectFactory as COF,
    Config,
    ConfigurableObject,
)  # E501: ignore
from jaix.runner import Runner, Optimiser
import logging
import copy
from gymnasium import Env
from typing import Type

logger = logging.getLogger("DefaultLogger")


class ATRunnerConfig(Config):
    def __init__(self, disp_interval: int = 20):
        self.disp_interval = disp_interval


class ATRunner(ConfigurableObject, Runner):
    config_class = ATRunnerConfig

    def run(self, env: Env, opt_class: Type[Optimiser], opt_config: Config):
        logger.debug("Starting experiment with %s on %s", opt_class, env)
        # Independent restarts (runs)
        while not env.stop():
            env.reset()
            logger.debug("Resetting optimiser")
            opt = COF.create(opt_class, opt_config, env=env)
            while not opt.stop() and not env.stop():
                X = opt.ask(env=env)
                res_list = []
                for x in X:
                    obs, r, term, trunc, info = env.step(x)
                    res_list.append(
                        {"obs": obs, "r": r, "term": term, "trunc": trunc, "info": info}
                    )
                # Reformat observations to dictlist
                # And pass as additional kwargs
                res_dict = {k: [dic[k] for dic in res_list] for k in res_list[0]}
                opt.tell(
                    env=env, solutions=X, function_values=res_dict["obs"], **res_dict
                )
                opt.disp(self.disp_interval)
                logger.debug(res_dict)

            info["opt_stop"] = opt.stop()
            info["env_stop"] = env.stop()
            logger.debug("Optimiser stopped.")
            logger.debug(
                f"Termination by opt {opt.stop()} env {env.stop()}"
            )  # TODO determine exact stopping criterion
            logger.debug(f"Result {info}")

        logger.debug("Experiment done")
