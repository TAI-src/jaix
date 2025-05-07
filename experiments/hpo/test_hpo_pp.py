from post_process import (
    bin_d,
    which2dec,
    which2bin,
    binary_distance,
    set_matches,
    mean_bin_dist,
    summarise,
)
import numpy as np


def test_bin_d():
    xi = [1, 0, 1, 0]
    xj = [0, 1, 1, 0]
    assert bin_d(xi, xj) == 2


def test_which2dec():
    x = [[2, 3]]
    assert which2dec(x) == [2**2 + 2**3]


def test_which2bin():
    x = [[2, 3]]
    n = 5
    assert which2bin(x, n) == [[0, 0, 1, 1, 0]]
    assert which2bin([[]], n) == [[0, 0, 0, 0, 0]]


def test_set_matches():
    xi = [[0, 1, 2], [1, 2, 3]]
    xj = [[0, 1], [1, 2]]
    assert set_matches(xi, xj) == 0
    xj = [[0, 1], [1, 2, 3]]
    assert set_matches(xi, xj) == 1


def test_binary_distance():
    xi = [[0, 1, 2], [1, 2, 3]]
    xj = [[0, 1], [1, 2]]
    assert binary_distance(xi, n=9, xj=xj) == [1, 1]
    xj = [[0, 1], [1, 2, 3]]
    assert binary_distance(xi, n=9, xj=xj) == [1, 0]
    xj = [[]]
    assert binary_distance(xi, n=9, xj=xj) == [3, 3]


def test_summarise():
    x = [1, 1, 1, 1]
    sum_dict = summarise(x)
    assert sum_dict["mean"] == 1
    assert sum_dict["min"] == 1
    assert sum_dict["max"] == 1
    assert sum_dict["std"] == 0
    assert sum_dict["med"] == 1

    sum_dict = summarise([])
    assert np.isnan(sum_dict["mean"])
    assert np.isnan(sum_dict["min"])
    assert np.isnan(sum_dict["max"])
    assert np.isnan(sum_dict["std"])
    assert np.isnan(sum_dict["med"])

    sum_dict = summarise([1])
    assert sum_dict["mean"] == 1
    assert sum_dict["min"] == 1
    assert sum_dict["max"] == 1
    assert sum_dict["std"] == 0
    assert sum_dict["med"] == 1


def test_mean_bin_dist():
    sum_dict = mean_bin_dist(1, 9)
    assert sum_dict["min"] == 0
    assert sum_dict["max"] == 9
    sum_dict = mean_bin_dist(2, 9)
    assert sum_dict["max"] == 8
