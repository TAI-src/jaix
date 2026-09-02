from jaix.utils.dirty_list import DirtyList


def test_dirty_list():
    dlist = DirtyList([1, 2, 3])
    assert dlist.dirty is True
    dlist.mark_clean()
    assert dlist.dirty is False
    dlist.append(4)
    assert dlist.dirty is True
    dlist.mark_clean()
    dlist.extend([5, 6])
    assert dlist.dirty is True

    dlist_empty = DirtyList([])
    assert dlist_empty.dirty is False
