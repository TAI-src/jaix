from ttex.config import ConfigurableObject, Config
import gymnasium as gym
import numpy as np
from jaix.env.singular import SingularEnvironment
from jaix.env.utils.ase import LJClustAdapter, LJClustAdapterConfig

from jaix import LOGGER_NAME
import logging

logger = logging.getLogger(LOGGER_NAME)

class LJClustEnvironmentConfig(Config):
    def __init__(self,
                 ljclust_adapter_config: LJClustAdapterConfig,
                 budget: int = np.iinfo(np.int32).max,
                 target_accuracy: float: 0.0,
                 ):
        self.ljclust_adapter_config = ljclust_adapter_config
        self.budget = budget
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

    def __init__(self, config: LJClustEnvironmentConfig,
                 func: int, inst:int):
        ConfigurableObject.__init__(self, config)
        SingularEnvironment.__init__(self, func, inst)
        species_str = LJClustAdapter.get_species_str(func, inst)
        self.adapter.set_species(species_str)

        # TODO need to figure out the actual box where to look for atom positions
        self.action_space = gym.spaces.Box(
            low=0.0,
            high=np.inf,
            shape=(self.adapter.num_atoms, 3),
            dtype=np.float64
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf,
            high=np.inf,
            shape=(1,),
            dtype=np.float64
        )

    def _get_info(self):
        return {
            "species": self.adapter.atom_str,
            "num_atoms": self.adapter.num_atoms,
            "box_length": self.adapter.box_length,
            "min_val": self.adapter.min_val,
        }

    def setup():
        pass
        
