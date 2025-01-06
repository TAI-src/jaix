from .. import DummyEnv
import numpy as np


def loop(dimension, num_objectives, opt, env=None):
    assert opt.stop() == {}
    if env is None:
        env = DummyEnv(dimension=dimension, num_objectives=num_objectives)
    # ask
    X = opt.ask(env=env)
    for x in X:
        assert env.action_space.contains(np.asarray(x, dtype=env.action_space.dtype))

    # tell
    Y = [env.step(x)[0] for x in X]  # First values is observation
    info = [{} for _ in range(len(Y))]
    r = [y[0] for y in Y]
    for y in Y:
        assert env.observation_space.contains(y)
    opt.tell(env=env, solutions=X, function_values=Y, info=info, r=r)

    return X, Y
