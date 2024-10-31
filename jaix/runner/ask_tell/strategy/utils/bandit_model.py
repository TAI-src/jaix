from ttex.config import Config, ConfigurableObject
import numpy as np
from enum import Enum

# https://medium.com/@enendufrankc/implement-a-multi-armed-bandit-algorithm-18afa2354c3f

BanditExploitStrategy = Enum("BanditExploitStrategy", ["MAX", "PROP"])


class BanditConfig(Config):
    def __init__(
        self,
        epsilon: float,
        min_tries: int = 10,
        exploit_strategy: BanditExploitStrategy = BanditExploitStrategy.MAX,
    ):
        self.epsilon = epsilon
        # Choices are tried at least min_tries times before exploiting
        self.min_tries = min_tries
        self.exploit_strategy = exploit_strategy


class Bandit(ConfigurableObject):
    config_class = BanditConfig

    def __init__(self, config: BanditConfig, num_choices: int):
        ConfigurableObject.__init__(self, config)
        """
        Initialises a multi-armed bandit problems with n arms.
        Q: The agent's estimated mean reward for each arm
        N: The number of times each arm has been chosen
        """
        self.n = num_choices
        self.Q = np.zeros(self.n)
        self.N = np.zeros(self.n)

    def next_choice(self):
        if np.random.rand() < self.epsilon:
            # Explore: Choose a random action
            choice = np.random.randint(len(self.Q))
        elif any(self.N <= self.min_tries):
            # Not all choices have been tried the minimum number of tries
            next_idx = np.where(self.N <= self.min_tries)[0]
            choice = np.random.choice(next_idx)
        elif self.exploit_strategy == BanditExploitStrategy.MAX:
            # Pick choice with highest Q
            next_idx = np.where(self.Q == np.max(self.Q))[0]
            choice = np.random.choice(next_idx)
        elif self.exploit_strategy == BanditExploitStrategy.PROP:
            # Pick choice with probability proportional to Q
            p = self.Q / sum(self.Q)
            choice = np.random.choice(range(len(self.Q)), p=p)
        else:
            raise ValueError(f"undefined {self.exploit_strategy}")
        return choice

    def update_stats(self, choice, reward):
        self.N[choice] += 1
        self.Q[choice] += (reward - self.Q[choice]) / self.N[choice]
