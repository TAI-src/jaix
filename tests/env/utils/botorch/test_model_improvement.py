from jaix.env.utils.botorch.model_improvement import (
    ModelImprovementEvaluator,
    ModelImprovementEvaluatorConfig,
    AcqFunc,
)
from jaix.env.utils.problem.static_problem import StaticProblem
import torch
from botorch.test_functions import Branin
from typing import List, Tuple


class BraninStaticProblem(StaticProblem):
    def __init__(self):
        self.lower_bounds = [-5.0, 0.0]
        self.upper_bounds = [10.0, 15.0]
        super().__init__(dimension=2, num_objectives=1)

    def __call__(self, x: List[float]) -> Tuple[List[float], List[float]]:
        # Branin function is defined for 2D input, so we need to ensure x is 2D
        assert len(x) == 2, "Input must be a 2D point"
        if isinstance(x, list):
            x = torch.tensor(x, dtype=torch.float32)
        else:
            raise ValueError("Input must be a list or a torch.Tensor")
        retval = Branin(negate=True)(x)
        # currently not doing noise, but we could add noise here if desired
        return [retval.item()], [
            retval.item()
        ]  # Return as lists to match expected output format


def get_config():
    return ModelImprovementEvaluatorConfig(
        acq_func=AcqFunc.MAX_VALUE_ENTROPY, num_initial_points=5
    )


def test_init_model():
    model = ModelImprovementEvaluator.init_model(
        BraninStaticProblem(), num_initial_points=5
    )
    assert model is not None
    assert model.train_inputs[0].shape == (5, 2)
    X = model.train_inputs[0]
    assert torch.all(X[:, 0] >= -5.0) and torch.all(X[:, 0] <= 10.0)
    assert torch.all(X[:, 1] >= 0.0) and torch.all(X[:, 1] <= 15.0)
    Y = model.train_targets
    model.eval()
    post = model.posterior(X)
    assert post.mean.shape == (5, 1)
    assert post.variance.shape == (5, 1)
    # compare the mean predictions to the actual function values at the training points
    assert torch.allclose(post.mean.squeeze(), Y.squeeze(), atol=1e-1)


def test_get_score_func():
    model = ModelImprovementEvaluator.init_model(
        BraninStaticProblem(), num_initial_points=5
    )
    candidate_set = torch.tensor(
        [[-5.0, 0.0], [10.0, 15.0], [2.5, 7.5]], dtype=torch.float32
    )
    score_func = ModelImprovementEvaluator.get_score_func(
        model, candidate_set, AcqFunc.MAX_VALUE_ENTROPY
    )
    scores = score_func([-5.0, 0.0])
    print(scores)
