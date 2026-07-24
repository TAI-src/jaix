from jaix.env.wrapper.value_track_wrapper import ValueTrackWrapper
import pytest
import numpy as np
import gymnasium as gym


def test_get_vals():
    obs = [5]
    r = 10
    info = {"val": 15}
    assert ValueTrackWrapper.get_vals(obs, r, info, "obs0") == {"obs0": 5}
    assert ValueTrackWrapper.get_vals(np.array(obs), r, info, ["obs0"]) == {"obs0": 5}
    assert ValueTrackWrapper.get_vals((0, np.array([obs])), r, info, ["obs0"]) == {
        "obs0": 5
    }
    assert ValueTrackWrapper.get_vals(obs, r, info, ["r"]) == {"r": 10}
    assert ValueTrackWrapper.get_vals(obs, r, info, ["val"]) == {"val": 15}
    with pytest.raises(AssertionError):
        ValueTrackWrapper.get_vals([5, 6], r, info, ["obs0"])
    with pytest.raises(ValueError):
        ValueTrackWrapper.get_vals(obs, r, info, ["invalid"])

    # Test with multiple state_eval keys
    assert ValueTrackWrapper.get_vals(obs, r, info, ["obs0", "r", "val"]) == {
        "obs0": 5,
        "r": 10,
        "val": 15,
    }


@pytest.mark.parametrize("state_eval", ["r", ["r", "r"]])
def test_val_updates(state_eval):
    env = gym.make("Pendulum-v1", render_mode="rgb_array")
    env.reset(seed=1337)
    wrapped_env = ValueTrackWrapper(env, state_eval=state_eval, is_min=True)

    wrapped_env.reset(seed=1337)
    assert wrapped_env.steps == 0
    assert len(wrapped_env.first_val) == 0
    assert len(wrapped_env.best_val) == 0
    assert len(wrapped_env.last_val) == 0

    env_cpy = gym.make("Pendulum-v1", render_mode="rgb_array")
    env_cpy.reset(seed=1337)

    num_steps = 1000
    rewards = []
    for i in range(num_steps):
        act = wrapped_env.action_space.sample()
        wrapped_env.step(act)
        assert wrapped_env.steps == i + 1
        obs, r, _, _, info = env_cpy.step(act)
        val = ValueTrackWrapper.get_vals(obs, r, info, state_eval)
        rewards.append(val)

    assert rewards[0] == wrapped_env.first_val
    for se in state_eval:
        # we assume min
        assert wrapped_env.best_val[se] == min([r[se] for r in rewards])
    assert rewards[-1] == wrapped_env.last_val
