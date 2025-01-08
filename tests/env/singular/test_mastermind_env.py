from jaix.env.singular import MastermindEnvironmentConfig, MastermindEnvironment
from gymnasium.utils.env_checker import check_env
import pytest
import copy


def test_init():
    config = MastermindEnvironmentConfig()
    env = MastermindEnvironment(config, inst=1)
    assert env.num_slots <= config.num_slots_range[1]
    assert env.num_slots >= config.num_colours_range[0]
    assert all(env.num_colours <= config.num_colours_range[1])
    assert all(env.num_colours >= config.num_colours_range[0])
    assert env.action_space.contains(env._solution)
    assert env.action_space.contains(env.num_colours - 1)
    assert not all(env.num_colours == env.num_colours[0])


def test_basic():
    config = MastermindEnvironmentConfig()
    env = MastermindEnvironment(config, inst=21)
    check_env(env)


@pytest.mark.parametrize("seq", [True, False])
def test_step_non_sequential(seq):
    config = MastermindEnvironmentConfig(sequential=seq, max_guesses=2)
    env = MastermindEnvironment(config, inst=3)
    obs, r, term, trunc, info = env.step(env._solution)
    assert obs[0] == r
    assert r == 0
    assert term
    assert not trunc

    all_wrong = env._solution + [1] * env.num_slots
    obs, r, term, trunc, info = env.step(all_wrong)
    assert r == env.num_slots
    assert not term
    assert trunc

    for i in range(env.num_slots):
        one_wrong = copy.deepcopy(env._solution)
        one_wrong[i] += 3
        obs, r, term, trunc, info = env.step(one_wrong)
        if seq:
            # Fitness depends on which one is wrong
            assert r == env.num_slots - i
        else:
            # Only one wrong, so fitness is 1
            assert r == 1
        assert not term
        assert trunc
