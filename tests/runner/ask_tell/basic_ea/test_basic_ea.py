from jaix.runner.ask_tell.strategy.basic_ea import BasicEA, BasicEAConfig
from jaix.runner.ask_tell.strategy.basic_ea import (
    EAStrategy,
    MutationOp,
    CrossoverOp,
    UpdateStrategy,
    WarmStartStrategy,
)
import pytest
from .. import DummyEnv, loop
import numpy as np
from gymnasium import spaces


@pytest.mark.parametrize(
    "strategy,mu,lam,mutation_op,crossover_op,mutation_opts,crossover_opts,low,high",
    [
        (
            EAStrategy.Comma,
            2,
            4,
            MutationOp.FLIP,
            CrossoverOp.ONEPOINT,
            {"p": 0.5},
            {"k": 0},
            2,
            5,
        ),
        (EAStrategy.Plus, 1, 1, MutationOp.FLIP, None, {}, {}, 2, 2),
        (EAStrategy.Plus, 10, 5, None, CrossoverOp.UNIFORM, {}, {}, 2, 2),
    ],
)
def test_basic(
    strategy,
    mu,
    lam,
    mutation_op,
    crossover_op,
    mutation_opts,
    crossover_opts,
    low,
    high,
):
    config = BasicEAConfig(
        strategy=strategy,
        mu=mu,
        lam=lam,
        mutation_op=mutation_op,
        crossover_op=crossover_op,
        mutation_opts=mutation_opts,
        crossover_opts=crossover_opts,
    )
    dimension = 3
    num_acts = np.random.randint(
        low=low,
        high=high + 1,
        size=dimension,
    )
    action_space = spaces.MultiDiscrete(num_acts)
    env = DummyEnv(num_objectives=1, action_space=action_space, dimension=dimension)
    ea = BasicEA(config, env)
    for _ in range(10):
        X = ea.ask(env)
        assert len(X) == lam
        for x in X:
            assert len(x) == dimension
            assert x in action_space
        ea.tell(env, X, [[np.random.rand()] for _ in range(lam)])
        assert len(ea.pop) == mu
        assert ea.pop[0].fitness <= ea.pop[-1].fitness


def test_assertions_config():
    # Test that lam <= mu for no crossover
    with pytest.raises(AssertionError):
        BasicEAConfig(
            strategy=EAStrategy.Comma,
            mu=2,
            lam=1,
            mutation_op=MutationOp.FLIP,
            crossover_op=None,
            mutation_opts={},
            crossover_opts={},
        )
    # Test that lam >= mu for Comma strategy
    with pytest.raises(AssertionError):
        BasicEAConfig(
            strategy=EAStrategy.Comma,
            mu=2,
            lam=1,
            mutation_op=MutationOp.FLIP,
            crossover_op=CrossoverOp.ONEPOINT,
            mutation_opts={},
            crossover_opts={},
        )
    # Test that at least one of mutation_op or crossover_op is not None
    with pytest.raises(AssertionError):
        BasicEAConfig(
            strategy=EAStrategy.Comma,
            mu=2,
            lam=3,
            mutation_op=None,
            crossover_op=None,
            mutation_opts={},
            crossover_opts={},
        )
    # Test that mu >= 2 for crossover
    with pytest.raises(AssertionError):
        BasicEAConfig(
            strategy=EAStrategy.Plus,
            mu=1,
            lam=1,
            mutation_op=None,
            crossover_op=CrossoverOp.UNIFORM,
            mutation_opts={},
            crossover_opts={},
        )


def test_improvement():
    config = BasicEAConfig(
        strategy=EAStrategy.Plus,
        mu=1,
        lam=1,
        mutation_op=MutationOp.FLIP,
        crossover_op=None,
        mutation_opts={},
        crossover_opts={},
    )
    dimension = 5
    action_space = spaces.MultiBinary(dimension)
    env = DummyEnv(dimension=dimension, action_space=action_space, num_objectives=1)
    ea = BasicEA(config, env)
    y_coll = []
    for i in range(10):
        if i == 0:
            init_fitness = min(ind.fitness for ind in ea.pop)
        _, Y = loop(3, 1, ea, env)
        y_coll.append(Y[0])
    final_fitness = min(ind.fitness for ind in ea.pop)
    assert final_fitness <= init_fitness
    assert final_fitness == min(y_coll)


def test_age():
    config = BasicEAConfig(
        strategy=EAStrategy.Comma,
        mu=1,
        lam=1,
        mutation_op=MutationOp.FLIP,
        crossover_op=None,
        mutation_opts={},
        crossover_opts={},
    )
    dimension = 5
    action_space = spaces.MultiBinary(dimension)
    env = DummyEnv(dimension=dimension, action_space=action_space, num_objectives=1)
    ea = BasicEA(config, env)
    for _ in range(10):
        loop(3, 1, ea, env)
    assert ea.pop[0].generation == 10


@pytest.mark.parametrize("warm_start_best", [True, False])
def test_warm_start(warm_start_best):
    config = BasicEAConfig(
        strategy=EAStrategy.Plus,
        mu=1,
        lam=1,
        mutation_op=MutationOp.FLIP,
        crossover_op=None,
        mutation_opts={},
        crossover_opts={},
    )
    if warm_start_best:
        config.warm_start_strategy = WarmStartStrategy.BEST
    else:
        config.warm_start_strategy = WarmStartStrategy.LAST
    dimension = 14
    action_space = spaces.MultiBinary(dimension)
    env = DummyEnv(dimension=dimension, action_space=action_space, num_objectives=1)
    ea = BasicEA(config, env)
    for _ in range(10):
        X, Y = loop(3, 1, ea, env)
    r1 = min(ind.fitness for ind in ea.pop)
    env.reset()
    xlast = X[-1]
    xfav = sorted(ea.pop, key=lambda x: x.fitness)[0].x
    ea.warm_start(xlast, env)
    if warm_start_best:
        assert all(ea.pop[0].x == xfav)
    else:
        assert all(ea.pop[0].x == xlast)
    for _ in range(10):
        X, Y = loop(3, 1, ea, env)
    r2 = min(ind.fitness for ind in ea.pop)
    if warm_start_best:
        assert r2 <= r1


def test_update():
    config = BasicEAConfig(
        strategy=EAStrategy.Plus,
        mu=1,
        lam=1,
        mutation_op=MutationOp.FLIP,
        crossover_op=None,
        mutation_opts={},
        crossover_opts={},
        update_strategy=UpdateStrategy.DDL,
        update_opts={"F": 1, "s": 1},
    )
    dimension = 5
    action_space = spaces.MultiBinary(dimension)
    env = DummyEnv(dimension=dimension, action_space=action_space, num_objectives=1)
    ea = BasicEA(config, env)
    assert "p" not in ea.mutation_opts
    for _ in range(10):
        loop(3, 1, ea, env)
    assert "p" in ea.mutation_opts
    assert ea.mutation_opts["p"] == 1 / dimension
