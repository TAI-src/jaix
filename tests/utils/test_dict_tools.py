from jaix.utils.dict_tools import nested_set


def test_nested_set():
    dict = {"a": {"b": {"c": 1, "d": 2}, "e": 3}, "f": 4}
    nested_set(dict, ["a", "b", "c"], 2)
    assert dict == {"a": {"b": {"c": 2, "d": 2}, "e": 3}, "f": 4}
