from jaix.env.singular.mastermind_env import (
    MastermindEnvironmentConfig,
    MastermindEnvironment,
)
from gymnasium.utils.env_checker import check_env
import pytest
import copy


def test_init():
    config = MastermindEnvironmentConfig()
    env = MastermindEnvironment(config, func=0, inst=1)
    assert env.num_slots <= config.num_slots_range[1]
    assert env.num_slots >= config.num_colours_range[0]
    assert all(env.num_colours <= config.num_colours_range[1])
    assert all(env.num_colours >= config.num_colours_range[0])
    assert env.action_space.contains(env._solution)
    assert env.action_space.contains(env.num_colours - 1)
    assert not all(env.num_colours == env.num_colours[0])


def test_basic():
    config = MastermindEnvironmentConfig()
    env = MastermindEnvironment(config, func=1, inst=21)
    check_env(env)


@pytest.mark.parametrize("seq", [True, False])
def test_step_non_sequential(seq):
    config = MastermindEnvironmentConfig(max_guesses=2)
    env = MastermindEnvironment(config, func=seq, inst=3)
    obs, r, term, trunc, info = env.step(env._solution)
    assert len(obs) == 1
    assert r is None
    assert term
    assert not trunc

    all_wrong = env._solution + [1] * env.num_slots
    obs, r, term, trunc, info = env.step(all_wrong)
    assert r is None
    assert obs[0] == env.num_slots
    assert not term
    assert trunc

    env._order = list(range(env.num_slots))

    for i in range(env.num_slots):
        one_wrong = copy.deepcopy(env._solution)
        one_wrong[i] += 3
        obs, r, term, trunc, info = env.step(one_wrong)
        assert len(obs) == 1
        if seq:
            # Fitness depends on which one is wrong
            assert obs[0] == env.num_slots - i
        else:
            # Only one wrong, so fitness is 1
            assert obs[0] == 1
        assert not term
        assert trunc


def test_order():
    config = MastermindEnvironmentConfig(max_guesses=2)
    env = MastermindEnvironment(config, func=True, inst=3)

    act = copy.deepcopy(env._solution)
    act[0] += 1

    env._order = list(range(env.num_slots))
    env._order.reverse()

    obs, _, _, _, _ = env.step(act)
    assert obs[0] == 1  # Because of reverse order


def test_inst_seeding():
    config = MastermindEnvironmentConfig()
    env1 = MastermindEnvironment(config, func=0, inst=1)
    env2 = MastermindEnvironment(config, func=0, inst=1)
    env3 = MastermindEnvironment(config, func=0, inst=2)

    assert env1.num_slots == env2.num_slots
    assert all(env2.num_colours == env2.num_colours)
    assert all(env1._solution == env2._solution)

    act = env1.action_space.sample()
    obs1, _, _, _, _ = env1.step(act)
    obs2, _, _, _, _ = env2.step(act)
    assert obs1 == obs2
    if env3.action_space.contains(act):
        obs3, _, _, _, _ = env3.step(act)
        assert obs1 != obs3
