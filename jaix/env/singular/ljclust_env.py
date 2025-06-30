from ttex.config import ConfigurableObject, Config, ConfigurableObjectFactory as COF
import gymnasium as gym
import numpy as np
from jaix.env.singular import SingularEnvironment
from jaix.env.utils.ase import LJClustAdapter, LJClustAdapterConfig

from jaix import LOGGER_NAME
import logging

logger = logging.getLogger(LOGGER_NAME)


class LJClustEnvironmentConfig(Config):
    def __init__(
        self,
        ljclust_adapter_config: LJClustAdapterConfig,
        target_accuracy: float = 1e-5,
    ):
        self.ljclust_adapter_config = ljclust_adapter_config
        self.target_accuracy = target_accuracy


class LJClustEnvironment(ConfigurableObject, SingularEnvironment):
    config_class = LJClustEnvironmentConfig

    @staticmethod
    def info(config: LJClustEnvironmentConfig):
        # Return information about the environment
        # TODO: Need to figure out what could be used for different functions and instances
        # then read from adapter instead of hardcoding
        return {
            "num_funcs": 1,
            "num_insts": 148,
        }

    def __init__(self, config: LJClustEnvironmentConfig, func: int, inst: int):
        ConfigurableObject.__init__(self, config)
        SingularEnvironment.__init__(self, func, inst)
        species_str = LJClustAdapter.finst2species(func, inst)
        self.adapter = COF.create(LJClustAdapter, config.ljclust_adapter_config)
        self.adapter.set_species(species_str)

        # TODO need to figure out the actual box where to look for atom positions
        self.action_space = gym.spaces.Box(
            low=0.0, high=np.inf, shape=(self.adapter.num_atoms, 3), dtype=np.float64
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf, high=np.inf, shape=(1,), dtype=np.float64
        )
        self.best_so_far = np.inf

    def _get_info(self):
        return {
            "species": self.adapter.atom_str,
            "num_atoms": self.adapter.num_atoms,
            "box_length": self.adapter.box_length,
            "min_val": self.adapter.min_val,
            "best_so_far": self.best_so_far,
        }

    def stop(self) -> bool:
        # Stop if the best energy is below the target accuracy
        return self.best_so_far - self.adapter.min_val <= self.target_accuracy

    def reset(self, seed=None, options=None):
        """
        Resets the environment to an initial state,
        required before calling step.
        Returns the first agent observation for an episode and information,
        i.e. metrics, debug info.
        """
        if options is None or "online" not in options or not options["online"]:
            # We only do partial resets for ec, so still "online"
            raise ValueError("HPO environments are always online")
        self.num_resets += 1
        self.best_so_far = np.inf
        return None, self._get_info()

    def step(self, pos):
        val: float
        if not self.adapter.validate(pos):
            val = np.inf  # TODO: Figure out proper penalties or nan
            add_info = {"invalid": True}
        else:
            val, add_info = self.adapter.evaluate(pos)
            add_info["invalid"] = False
            if val < self.best_so_far:
                self.best_so_far = val
        info = self._get_info()
        info.update(add_info)
        logger.debug(f"Step: {pos}, Info: {info}")
        terminated = self.stop()
        truncated = False
        r = -val
        return val, r, terminated, truncated, info

    def render(self):
        """
        Render the environment.
        This method is not implemented as rendering is not required for this environment.
        """
        logger.debug(self._get_info())
