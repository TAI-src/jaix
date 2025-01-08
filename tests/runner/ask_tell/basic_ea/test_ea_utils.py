from jaix.runner.ask_tell.strategy.utils.ea_utils import (
    global_flip,
    onepoint_crossover,
    uniform_crossover,
    Individual,
    select,
)
import pytest
import numpy as np


def test_global_flip_deterministic():
    # test binary deterministic
    x = [1, 0, 0, 1, 0]
    p = 1
    x = global_flip(x, p)
    assert x == [0, 1, 1, 0, 1]

    # test integer deterministic
    x = [1, 2, -3, 4, 5]
    xn = global_flip(x, p, low=-3, high=5)
    assert all([vx != vn for vx, vn in zip(x, xn)])


@pytest.mark.parametrize(
    "p, low, high",
    [
        (0.5, 0, 1),
        (0.5, 1, 5),
    ],
)
def test_global_flip_stochastic(p, low, high):
    # test binary stochastic
    x = np.random.randint(low, high + 1, 5)
    diffs = 0
    reps = 100
    for _ in range(reps):
        xn = global_flip(x, p, low, high)
        diffs += sum([vx != vn for vx, vn in zip(x, xn)])
    assert diffs < (p + 0.05) * reps * len(x)
    assert diffs > (p - 0.05) * reps * len(x)


def test_global_flip_diff_values():
    # deterministic
    x = [5, 4, 3, 2, 1]
    low = [4, 3, 2, 1, 0]
    high = [5, 4, 3, 2, 1]
    p = 1
    xn = global_flip(x, p, low=low, high=high)
    assert all([vx != vn for vx, vn in zip(x, xn)])
    assert all(np.array(xn) == low)

    xn = global_flip(x, p, low=low, high=5)
    assert xn[0] == 4
    assert not all(xn == xn[0])


def test_onepoint_crossover():
    # test deterministic
    x1 = [1, 2, 3, 4, 5]
    x2 = [6, 7, 8, 9, 10]
    xc = onepoint_crossover(x1, x2, k=2)
    assert all(xc == [1, 2, 8, 9, 10])
    with pytest.raises(AssertionError):
        xc = onepoint_crossover(x1, x2, k=5)
    xc = onepoint_crossover(x1, x2, k=0)
    assert all(xc == [6, 7, 8, 9, 10])

    # test stochastic
    xc = onepoint_crossover(x1, x2)
    for i in range(len(x1)):
        assert xc[i] in [x1[i], x2[i]]

    # test different lengths
    x1 = [1, 2, 3]
    x2 = [4, 5]
    with pytest.raises(AssertionError):
        xc = onepoint_crossover(x1, x2)


def test_uniform_crossover():
    # test deterministic
    x1 = [1, 2, 3, 4, 5]
    x2 = [6, 7, 8, 9, 10]
    mask = [1, 0, 1, 0, 1]
    xc = uniform_crossover(x1, x2, mask=mask)
    assert all(xc == [1, 7, 3, 9, 5])

    # test stochastic
    xc = uniform_crossover(x1, x2)
    for i in range(len(x1)):
        assert xc[i] in [x1[i], x2[i]]

    # test different lengths
    x1 = [1, 2, 3]
    x2 = [4, 5]
    with pytest.raises(AssertionError):
        xc = uniform_crossover(x1, x2)

    # test invalid mask
    mask = [1, 0, 2, 0, 1]
    with pytest.raises(AssertionError):
        xc = uniform_crossover(x1, x2, mask=mask)
    mask = [1, 0, 1, 0]
    with pytest.raises(AssertionError):
        xc = uniform_crossover(x1, x2, mask=mask)


def test_select():
    # test deterministic
    pop = [
        Individual([1, 2, 3], 1, 0),
        Individual([4, 5, 6], 2, 0),
        Individual([7, 8, 9], 3, 0),
    ]
    pop.reverse()
    mu = 2
    selected = select(pop, mu)
    assert len(selected) == mu
    assert selected[0].fitness == 1
    assert selected[1].fitness == 2
    assert selected[0].x == [1, 2, 3]
    assert selected[1].x == [4, 5, 6]
