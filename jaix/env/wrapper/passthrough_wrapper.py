import gymnasium as gym
from typing import Dict


class PassthroughWrapper(gym.Wrapper):
    def __init__(self, env: gym.Env, passthrough: bool):
        super().__init__(env)
        self.passthrough = passthrough

    def __getattr__(self, name):
        # Pass through anything else that is not overriden
        # This function is only called after it is not found in self
        if self.passthrough:  # Passing through
            return getattr(self.env, name)
        else:  # Stay consistent with default behaviour
            raise AttributeError(name)

    def _stop(self) -> Dict:
        return {}

    def _rec_stop(self) -> Dict:
        term_conds = self._stop()
        if isinstance(self.env, PassthroughWrapper):
            term_conds.update(self.env._rec_stop())
        return term_conds

    def stop(self) -> Dict:
        term_conds = self._rec_stop()
        if hasattr(self.env.unwrapped, "stop") and self.env.unwrapped.stop():
            term_conds["unwrapped"] = self.env.unwrapped.stop()
        return term_conds
