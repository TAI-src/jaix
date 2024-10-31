from ttex.config import (
    Config,
    ConfigurableObject,
    ConfigurableObjectFactory as COF,
)
from typing import List, Dict
from jaix.runner.ask_tell import ATOptimiser, ATOptimiserConfig, ATStrategy
from jaix.runner.ask_tell.strategy.utils import BanditConfig, Bandit
from gymnasium import Env


class ATBanditConfig(Config):
    def __init__(
        self,
        opt_confs: List[ATOptimiserConfig],
        bandit_config: BanditConfig,
    ):
        self.opt_confs = opt_confs
        self.bandit_config = bandit_config
        self.init_pop_size = opt_confs[0].init_pop_size


class ATBandit(ConfigurableObject, ATStrategy):
    config_class = ATBanditConfig

    def __init__(self, config: ATBanditConfig, env):
        ConfigurableObject.__init__(self, config)
        ATStrategy.__init__(self, None)
        # Initialise bandit
        num_choices = len(self.opt_confs)
        self.bandit = COF.create(Bandit, self.bandit_config, num_choices)

        # Initialise the first optimiser
        self.opt = COF.create(ATOptimiser, self.opt_confs[0], env)
        self._active_opt = 0

    def initialize(self):
        pass

    def _run_bandit(self, rewards):
        for r in rewards:
            self.bandit.update_stats(self._active_opt, r)
        if len(rewards) > 0:
            # There was a q update, choose new opt
            self._active_opt = self.bandit.next_choice()
        # Return if there was an update
        return len(rewards) > 0

    def _update(self, infos, env):
        final_rewards = [info["final_r"] for info in infos if "final_r" in info]
        updated = self._run_bandit(final_rewards)
        if updated:
            # Need to init new algorithm
            self.opt = COF.create(
                ATOptimiser,
                self.opt_confs[self._active_opt],
                env,
            )
        return updated

    def ask(self, env, **optional_kwargs):
        if self.opt.stop():
            # Algorithm wants to stop, reset the env
            # and update internal bandit
            _, info = env.reset(options={"online": True})
            self._update([info], env)
        # Ask the active optimiser
        return self.opt.ask(env, **optional_kwargs)

    def tell(
        self,
        solutions,
        function_values,
        info: Dict,
        r: float,
        env: Env,
        **optional_kwargs,
    ):
        # Update Q info if any env stopped
        updated = self._update(info, env)
        # Otherwise we would start with tell, not ask
        # And CMA-ES does not like that apparently
        # TODO: Figure out a good strategy to not lose info
        # For now throwing it away because otherwise might not be enough info
        # also not independent otherwise
        if not updated:
            return self.opt.tell(
                env=env,
                solutions=solutions,
                function_values=function_values,
                info=info,
                r=r,
                **optional_kwargs,
            )

    def stop(self):
        return False