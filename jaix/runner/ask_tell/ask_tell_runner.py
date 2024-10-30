from ttex.config import (
    ConfigurableObjectFactory as COF,
    Config,
    ConfigurableObject,
)  # E501: ignore
from jaix.runner import Runner, Optimiser
import logging
import gymnasium as gym
from typing import Type
from jaix.env.wrapper import (
    MaxEvalWrapper,
    MaxEvalWrapperConfig,
    WrappedEnvFactory as WEF,
)

logger = logging.getLogger("DefaultLogger")


class ATRunnerConfig(Config):
    def __init__(self, max_evals: int, disp_interval: int = 20):
        self.disp_interval = disp_interval
        self.max_evals = max_evals


class ATRunner(ConfigurableObject, Runner):
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
        wrappers = [(MaxEvalWrapper, MaxEvalWrapperConfig(max_evals=self.max_evals))]  # type: List[Tuple[Type[gym.Wrapper], Union[Config, Dict]]]

        wenv = WEF.wrap(env, wrappers)  # type: PassthroughWrapper
        # Independent restarts (runs)
        while not wenv.stop():
            wenv.reset()
            logger.debug("Resetting optimiser")
            opt = COF.create(opt_class, opt_config, env=wenv)
            while not opt.stop() and not wenv.stop():
                X = opt.ask(env=wenv)
                res_list = []
                for x in X:
                    obs, r, term, trunc, info = wenv.step(x)
                    res_list.append(
                        {"obs": obs, "r": r, "term": term, "trunc": trunc, "info": info}
                    )
                # Reformat observations to dictlist
                # And pass as additional kwargs
                res_dict = {k: [dic[k] for dic in res_list] for k in res_list[0]}
                opt.tell(
                    env=wenv, solutions=X, function_values=res_dict["obs"], **res_dict
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

        logger.debug("Experiment done")
