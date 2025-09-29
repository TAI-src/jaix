from jaix.runner.ask_tell.strategy.random import RandomATStrat, RandomATStratConfig
from ttex.config import ConfigFactory as COF
from jaix.runner.ask_tell.strategy.bandit import ATBandit, ATBanditConfig
from jaix.runner.ask_tell.at_optimiser import ATOptimiserConfig
from jaix.runner.ask_tell.strategy.bandit import BanditConfig
from jaix.runner.ask_tell.strategy.utils.bandit_model import BanditExploitStrategy
from jaix.env.wrapper.auto_reset_wrapper import AutoResetWrapper, AutoResetWrapperConfig
from . import DummyEnv, loop
import numpy as np


def get_bandit(num_choices: int = 2, stop_after: int = -1):
    opt_confs = []
    for _ in range(num_choices):
        config = ATOptimiserConfig(
            strategy_class=RandomATStrat,
            strategy_config=RandomATStratConfig(ask_size=5),
            init_pop_size=1,
        )
        if stop_after > 0:
            config.stop_after = stop_after
        opt_confs.append(config)

    bandit_config = BanditConfig(
        epsilon=0.1,
        min_tries=10,
        exploit_strategy=BanditExploitStrategy.MAX,
    )
    config = ATBanditConfig(
        opt_confs=opt_confs,
        bandit_config=bandit_config,
    )
    bandit = ATBandit(config, DummyEnv())
    return bandit


def test_init():
    bandit = get_bandit()
    assert isinstance(bandit.opt.strategy, RandomATStrat)
    assert bandit._active_opt == 0


def test_update():
    env = DummyEnv()
    bandit_opt = get_bandit(num_choices=3)
    _, info = env.reset(options={"online": True})

    # Checking init state
    assert sum(bandit_opt.bandit.N) == 0
    assert bandit_opt._active_opt == 0
    assert bandit_opt.opt.countiter == 0
    current_opt = bandit_opt.opt

    # No final r
    updated = bandit_opt._update([info] * 3, env)
    assert not updated
    assert bandit_opt.opt.countiter == 0
    assert bandit_opt.opt == current_opt
    assert sum(bandit_opt.bandit.N) == 0

    # Final r
    info["final_r"] = 1
    updated = bandit_opt._update([info] * 3, env)
    assert updated
    assert bandit_opt.opt.countiter == 0
    assert bandit_opt.opt != current_opt
    assert sum(bandit_opt.bandit.N) == 3
    assert bandit_opt.bandit.N[0] == 3


def test_tell_update():
    env = DummyEnv()
    bandit_opt = get_bandit(num_choices=2)

    # Checking init state
    assert sum(bandit_opt.bandit.N) == 0
    assert bandit_opt._active_opt == 0
    assert bandit_opt.opt.countiter == 0

    obs, info = env.reset(options={"online": True})

    x = bandit_opt.ask(env)
    y = [env.step(x) for x in x]
    r = 1

    assert sum(bandit_opt.bandit.N) == 0

    # Tell without final r - just regular execution
    bandit_opt.tell(
        x,
        y,
        obs=env.observation_space.sample(),
        env=env,
        r=[r],
        trunc=[False],
        term=[False],
        info=[info],
    )
    assert sum(bandit_opt.bandit.N) == 0
    assert bandit_opt.opt.countiter == 1

    info["final_obs"] = env.observation_space.sample()
    bandit_opt.tell(
        x,
        y,
        obs=env.observation_space.sample(),
        env=env,
        r=y,
        trunc=[False],
        term=[False],
        info=[info],
    )
    assert sum(bandit_opt.bandit.N) == 0
    assert bandit_opt.opt.countiter == 2

    info["final_r"] = r
    bandit_opt.tell(
        x,
        y,
        obs=env.observation_space.sample(),
        env=env,
        r=y,
        trunc=[False],
        term=[False],
        info=[info],
    )
    assert bandit_opt.bandit.N[0] == 1
    assert bandit_opt.bandit.Q[0] == 1
    assert bandit_opt.opt.countiter == 0

    info.pop("final_r")
    bandit_opt.tell(
        x,
        y,
        obs=env.observation_space.sample(),
        env=env,
        r=y,
        trunc=[False],
        term=[False],
        info=[info],
    )
    assert bandit_opt.bandit.N[0] == 1
    assert bandit_opt.bandit.Q[0] == 1
    assert bandit_opt.opt.countiter == 1


def test_ask_update():
    env = DummyEnv()
    bandit_opt = get_bandit(num_choices=2, stop_after=2)
    env = AutoResetWrapper(AutoResetWrapperConfig(), DummyEnv())

    obs, info = env.reset(options={"online": True})

    # The inner optimiser is supposed to stop after 2 tells
    for _ in range(2):
        bandit_opt.ask(env)
        obs, r, trunc, term, info = env.step([1, 1, 1])
        assert not trunc
        assert not term
        assert not env.stop()
        bandit_opt.tell(
            [[1, 1, 1]],
            [obs],
            obs=[obs],
            env=env,
            r=[r],
            trunc=[trunc],
            term=[term],
            info=[info],
        )
    assert bandit_opt.opt.stop()
    assert not env.unwrapped.stop()
    assert bandit_opt.bandit.N[0] == 0
    assert bandit_opt.opt.countiter == 2
    assert env.man_resets == 1

    # Assure that if opt decided to stop, everything is reset
    bandit_opt.ask(env)
    assert not bandit_opt.stop()
    assert not bandit_opt.opt.stop()
    assert bandit_opt.bandit.N[0] == 1
    assert bandit_opt.opt.countiter == 0
    assert env.man_resets == 2


def test_warm_start():
    env = DummyEnv()
    bandit_opt = get_bandit(num_choices=2, stop_after=2)
    env = AutoResetWrapper(AutoResetWrapperConfig(), DummyEnv())
    for _ in range(10):
        X, Y = loop(2, 1, bandit_opt, env)
    updated = bandit_opt.warm_start(None, env, None)
    assert bandit_opt._prev_r == [np.nan]
    assert updated
