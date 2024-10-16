from jaix.env import ResetStrategy
import gymnasium as gym
import logging

logger = logging.getLogger("DefaultLogger")
# TODO: Also use this wrapper to control randomness
# TODO: And to force reset before use

class JaixWrapper(gym.Wrapper):
    def __init__(self, env: gym.Env, default_reset_strat: ResetStrategy):
        self.default_reset_strat = default_reset_strat
        # Count how often the environment is reset
        # (corresponds to algorithms restarts +1 )
        self.num_resets = 0
        self.closed = False
        
    def reset(self, **kwargs):
        obs, info = self.env.reset(**kwargs)
        if self.steps >= self.min_steps:
            # Only update final r if it is not a failed reset
            info["final_r"] = self.prev_r
        self.man_resets += 1
        self.steps = 0
        self.prev_r = None
        return obs, info

    @property
    def state(self):
        """Simple representation of environment state
           To be implemented in specific environment
        """
        state = {"num_resets": self.num_resets, "close": self.closed}
        if hasattr(env, "state"):
            state.update(env.state)
        return state


    def render(self):
        """
        Renders the environments to help visualise what the agent see,
        examples modes are “human”, “rgb_array”, “ansi” for text.
        """
        logger.debug(self.state)



class ResetStrategy(Enum):
    BREAK = 0
    CONTINUE = 1


            if self.steps == 0:
                # This means we reset previously and on the first step
                # we are done. That is a fail
                self.failed_resets += 1
            elif self.steps >= self.min_steps:
                # Only update final r if it is not a failed reset
                info["final_r"] = r
            self.steps = 0
            self.auto_resets += 1
        else:
            self.steps += 1
        self.prev_r = r
        return obs, r, term, trunc, info

    def stop(self, failed_resets_thresh: int = 1):
        # Consider finally stopped if more failed resets than desired
        return self.failed_resets >= failed_resets_thresh
