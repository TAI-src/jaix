from ttex.config import ConfigurableObject, Config, ConfigurableObjectFactory as COF
import gymnasium as gym
import numpy as np
from jaix.env.singular.singular_environment import SingularEnvironment
from jaix.env.utils.ase import LJClustAdapter, LJClustAdapterConfig
from typing import Optional

import jaix.utils.globals as globals

import logging

logger = logging.getLogger(globals.LOGGER_NAME)


class LJClustEnvironmentConfig(Config):
    def __init__(
        self,
        ljclust_adapter_config: LJClustAdapterConfig,
        target_accuracy: float = 1e-5,
        by_species: bool = True,  # If True, function is species index, instance is number of atoms; otherwise, vice versa.
    ):
        self.ljclust_adapter_config = ljclust_adapter_config
        self.target_accuracy = target_accuracy
        self.by_species = by_species


class LJClustEnvironment(ConfigurableObject, SingularEnvironment):
    config_class = LJClustEnvironmentConfig

    @staticmethod
    def info(config: LJClustEnvironmentConfig):
        # Return information about the environment
        # TODO: Need to figure out what could be used for different functions and instances
        info = LJClustAdapter.get_info(config.by_species)
        return {
            "num_funcs": info["num_funcs"],
            "num_insts": info["num_insts"],
        }

    def __init__(self, config: LJClustEnvironmentConfig, func: int, inst: int):
        ConfigurableObject.__init__(self, config)
        SingularEnvironment.__init__(self, func, inst)
        species_str = LJClustAdapter.finst2species(func, inst)
        self.adapter = COF.create(LJClustAdapter, config.ljclust_adapter_config)
        self.adapter.set_species(species_str)

        # TODO: need to figure out the actual box where to look for atom positions
        # Based on adapter box length
        # An action is the positions of atoms in 3D space,
        # one coordinate per atom
        self.action_space = gym.spaces.Box(
            low=-np.inf,
            high=np.inf,
            shape=(self.adapter.num_atoms * 3,),
            dtype=np.float64,
        )
        # An observation is the energy after the last action
        self.observation_space = gym.spaces.Box(
            low=-np.inf, high=np.inf, shape=(1,), dtype=np.float64
        )
        self.best_so_far = np.inf
        # Count how often the environment is reset
        # (corresponds to algorithms restarts +1 )
        self.num_resets = 0

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
            raise ValueError("LJClustEnvironments are always online")
        self.num_resets += 1
        self.best_so_far = np.inf
        return None, self._get_info()

    def step(self, pos):
        val: Optional[float]
        pos = np.reshape(pos, (self.adapter.num_atoms, 3))
        if not self.adapter.validate(pos):
            val = None
            add_info = {"invalid": True}
        else:
            val, add_info = self.adapter.evaluate(pos)
            assert val is not None
            add_info["invalid"] = False
            if val < self.best_so_far:
                self.best_so_far = val
        info = self._get_info()
        info.update(add_info)
        logger.debug(f"Step: {pos}, Info: {info}")
        terminated = self.stop()
        truncated = False
        return (
            np.asarray([val], dtype=self.action_space.dtype),
            None,
            terminated,
            truncated,
            info,
        )

    def render(self):
        """
        Render the environment.
        This method is not implemented as rendering is not required for this environment.
        """
        logger.debug(self._get_info())
