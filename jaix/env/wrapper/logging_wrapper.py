import gymnasium as gym
import logging

logger = logging.getLogger("DefaultLogger")

class LoggingWrapper(gym.Wrapper):

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

