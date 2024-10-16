"""An Experiment"""
from ttex.config import (
    ConfigurableObjectFactory as COF, Config, ConfigurableObject
)  # E501: ignore
from jaix.runner import Runner
import logging
import copy

logger = logging.getLogger("DefaultLogger")


class ATRunnerConfig(Config):
    def __init__(self):
    	self.disp_interval = disp_interval

class ATRunner(ConfigurableObject, Runner):
    config_class = ATRunnerConfig

    def run(self, env, opt_class, opt_config):
        logger.debug("Starting experiment with %s on %s", opt_class, env)
        # Independent restarts (runs)
        while not env.stop():
            # BlackBox setting only makes sense with online resets
            env.reset(options={"online": True})
            logger.debug("Resetting optimiser")
            init_pop = env.unwrapped.sample_pop(opt_config.init_pop_size)
            opt = COF.create(opt_class, opt_config, init_pop)
            while not opt.stop() and not env.unwrapped.stop():
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
                # TODO: rewards should be single objective
                # Need to decide what values bb opt should get
                # Should define an action and observation space for algorithms
                # As well as potentially a checker for which wrappers
                opt.tell(X, res_dict["obs"], **res_dict, env=env)
                opt.disp(self.disp_interval)
                logger.debug(res_dict)


            info["opt_stop"] = opt.stop()
            info["env_stop"] = env.stop()
            logger.debug("Optimiser stopped.")
            logger.debug(
                f"Termination by opt {opt.stop()} env {env.unwrapped.stop()}"
            )  # TODO determine exact stopping criterion
            logger.debug(
                f"Result {info}"
            )  # TODO should not be calling private function

        logger.debug("Experiment done")
