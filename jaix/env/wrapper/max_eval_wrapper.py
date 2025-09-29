import gymnasium as gym
from ttex.config import ConfigurableObject, Config
from jaix.env.wrapper.passthrough_wrapper import PassthroughWrapper
from typing import Optional


class MaxEvalWrapperConfig(Config):
    def __init__(self, max_evals: int, passthrough: bool = True):
        self.max_evals = max_evals
        self.passthrough = passthrough


class MaxEvalWrapper(PassthroughWrapper, ConfigurableObject):
    config_class = MaxEvalWrapperConfig

    def __init__(self, config: MaxEvalWrapperConfig, env: gym.Env):
        ConfigurableObject.__init__(self, config)
        PassthroughWrapper.__init__(self, env, self.passthrough)
        self._evals = 0

    def reset(
        self,
        *,
        seed: Optional[int] = None,
        options: Optional[dict] = None,
    ):
        if options is None or "online" not in options or not options["online"]:
            # Reset evaluations _unless_ it is online
            self._evals = 0
        return self.env.reset(seed=seed, options=options)

    def step(self, action):
        self._evals += 1
        return self.env.step(action)

    def _stop(self):
        if self._evals >= self.max_evals:
            return {"max_evals": self._evals}
        else:
            return {}
