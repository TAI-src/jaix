from jaix.env.wrapper.online_wrapper import OnlineWrapper
from ttex.config import (
    ConfigurableObjectFactory as COF,
    Config,
)  # E501: ignore
from jaix.runner.runner import Runner
from jaix.runner.optimiser import Optimiser
import logging
import gymnasium as gym
from typing import Type, List, Tuple, Union, Dict
from jaix.env.wrapper.max_eval_wrapper import MaxEvalWrapper, MaxEvalWrapperConfig
from jaix.env.wrapper.passthrough_wrapper import PassthroughWrapper
from jaix.env.wrapper.wrapped_env_factory import WrappedEnvFactory as WEF

from jaix.utils.globals import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)


class ATRunnerConfig(Config):
    def __init__(self, max_evals: int, disp_interval: int = 20):
        self.disp_interval = disp_interval
        self.max_evals = max_evals


class ATRunner(Runner):
    config_class = ATRunnerConfig

    def run(
        self,
        env: gym.Env,
        opt_class: Type[Optimiser],
        opt_config: Config,
        *args,
        **kwargs,
    ):
        logger.debug("Starting experiment with %s on %s", opt_class, env)
        wrappers = [
            (MaxEvalWrapper, MaxEvalWrapperConfig(max_evals=self.max_evals)),
            (OnlineWrapper, {"online": True}),
        ]  # type: List[Tuple[Type[gym.Wrapper], Union[Config, Dict]]]

        wenv = WEF.wrap(env, wrappers)  # type: PassthroughWrapper
        # Independent restarts (runs)
        wenv.reset()
        while not wenv.stop():
            prev_id = wenv.id
            logger.debug("Resetting optimiser")
            opt = COF.create(opt_class, opt_config, env=wenv)
            logger.debug("Optimiser created")
            info = {}  # type: Dict
            while not opt.stop() and not wenv.stop():
                X = opt.ask(env=wenv)
                res_list = []
                for x in X:
                    logger.debug(f"Optimising {x}")
                    if wenv.id == prev_id:
                        # If the environment switches, the optimiser is reset
                        obs, r, term, trunc, info = wenv.step(x)
                        res_list.append(
                            {
                                "obs": obs,
                                "r": r,
                                "term": term,
                                "trunc": trunc,
                                "info": info,
                            }
                        )
                    else:
                        # Continue outer loop so that optimiser is reset to adjust to new action space
                        logger.debug(
                            f"Environment changed, warm_start optimiser for action_space {wenv.action_space} and env {wenv.id}"
                        )
                        prev_id = wenv.id
                        opt.warm_start(env=wenv, xstart=x, res_list=res_list)
                        break
                else:
                    # Reformat observations to dictlist
                    # And pass as additional kwargs
                    logger.debug(f"Res list {res_list}")
                    res_dict = {k: [dic[k] for dic in res_list] for k in res_list[0]}
                    opt.tell(
                        env=wenv,
                        solutions=X,
                        function_values=res_dict["obs"],
                        **res_dict,
                    )
                    opt.disp(self.disp_interval)
                    logger.debug(res_dict)

            info["opt_stop"] = opt.stop()
            info["env_stop"] = wenv.stop()
            logger.debug("Optimiser stopped.")
            logger.debug(
                f"Termination by opt {opt.stop()} env {wenv.stop()}"
            )  # TODO determine exact stopping criterion
            logger.debug(f"Result {info}")
            wenv.reset()

        logger.debug("Experiment done")
