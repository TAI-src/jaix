from jaix.suite.coco import COCOProblem
import pytest

if COCOProblem is not None:
    import cocoex as ex


@pytest.fixture(scope="session", autouse=True)
def skip_remaining_tests():
    try:
        import cocoex  # noqa: F401

        assert COCOProblem is not None
    except ImportError:
        assert COCOProblem is None
        pytest.skip(
            "Skipping COCO tests. If this is unexpected, check that the coco extra is installed."
        )


@pytest.fixture
def coco_problem():
    suite = ex.Suite("bbob-constrained", "", "")
    f = suite.get_problem(33)
    coco_problem = COCOProblem(f)
    return coco_problem


def test_coco_problem_init(coco_problem):
    # check attributes and init
    assert coco_problem.name is not None
    assert "bbob-constrained" in coco_problem.name
    init_sol = coco_problem.initial_solution_proposal()
    assert init_sol is not None


def test_evals_left(coco_problem):
    assert coco_problem.dimension == 2
    assert coco_problem.evalsleft(2) == 4
    coco_problem([0, 0])
    assert coco_problem.evalsleft(2) == 3
    coco_problem.constraint([0, 0])
    assert coco_problem.evalsleft(2) == 3
    coco_problem.constraint([0, 0])
    assert coco_problem.evalsleft(2) == 2


@pytest.mark.parametrize("n_obj", [1, 2])
def test_format(n_obj):
    if n_obj == 1:
        suite = ex.Suite("bbob", "", "dimensions:3")
    elif n_obj == 2:
        suite = ex.Suite("bbob-biobj", "", "dimensions:3")
    f = suite.get_problem(1)
    coco_problem = COCOProblem(f)
    res, _ = coco_problem([1, 2, 3])
    assert len(res) == n_obj
