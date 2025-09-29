import numpy as np
from jaix.runner.ask_tell.strategy.utils.bandit_model import (
    BanditConfig,
    Bandit,
    BanditExploitStrategy,
)
import pytest


def get_bandit(
    num_choices: int = 2,
    epsilon: float = 0.1,
    min_tries: int = 10,
    exploit_strategy: BanditExploitStrategy = BanditExploitStrategy.MAX,
):
    config = BanditConfig(
        epsilon=epsilon,
        min_tries=min_tries,
        exploit_strategy=exploit_strategy,
    )
    bandit = Bandit(config, num_choices=num_choices)
    return bandit


def test_init():
    bandit = get_bandit(num_choices=7, exploit_strategy=BanditExploitStrategy.PROP)
    assert len(bandit.Q) == 7


@pytest.mark.parametrize(
    "strategy", [BanditExploitStrategy.MAX, BanditExploitStrategy.PROP]
)
def test_bandit(strategy):
    bandit = get_bandit(num_choices=3, exploit_strategy=strategy)

    qstar = [1, 2, 3]
    num_steps = 1000
    successes = 0
    for _ in range(num_steps):
        action = bandit.next_choice()
        r = np.random.normal(qstar[action], 1)
        bandit.update_stats(action, r)
        successes += action == np.argmax(qstar)
    # Bandit should be able to identify the best option
    assert np.argmax(qstar) == np.argmax(bandit.Q)
    # Bandit should be better than random
    assert successes / num_steps >= 1 / len(qstar)
