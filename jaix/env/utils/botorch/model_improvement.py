from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.acquisition.max_value_entropy_search import qMaxValueEntropy
from ttex.config import Config, ConfigurableObject
from typing import List, Callable, Optional
import torch
from enum import Enum
from botorch.utils.transforms import normalize, standardize
from jaix.env.utils.problem.static_problem import StaticProblem
from botorch.acquisition import UpperConfidenceBound
from botorch.optim import optimize_acqf


class AcqFunc(Enum):
    MAX_VALUE_ENTROPY = "max_value_entropy"


class ModelImprovementEvaluatorConfig(Config):
    def __init__(
        self,
        acq_func: AcqFunc,
        num_initial_points: int,
        acq_func_params: dict | None = None,
        num_uniform_candidates: int = 100,
        num_ucb_candidates: int = 20,
        ucb_params: dict | None = None,
    ):
        Config.__init__(self)
        self.acq_func = acq_func
        self.acq_func_params = acq_func_params or {"num_fantasies": 20}
        self.num_initial_points = num_initial_points
        self.num_uniform_candidates = num_uniform_candidates
        self.num_ucb_candidates = num_ucb_candidates
        self.ucb_params = ucb_params or {
            "beta": 2.0,
            "num_restarts": 5,
            "raw_samples": 20,
        }


class ModelImprovementEvaluator(ConfigurableObject):
    config_class = ModelImprovementEvaluatorConfig

    def __init__(
        self,
        config: ModelImprovementEvaluatorConfig,
        func: StaticProblem,
        seed: int,
    ):
        ConfigurableObject.__init__(self, config)
        self.func = func
        self.seed = seed
        torch.manual_seed(seed)
        self.model = ModelImprovementEvaluator.init_model(
            func, self.config.num_initial_points
        )

    @staticmethod
    def init_model(func: StaticProblem, num_initial_points: int) -> SingleTaskGP:
        """
        Initialize a Gaussian Process model with initial data sampled from the function.
        The initial data is normalized and standardized before fitting the model.
        Args:
            func: The static problem function to sample from.
            num_initial_points: The number of initial data points to sample for model fitting.
        Returns:
            A fitted SingleTaskGP model based on the initial sampled data.
        """
        # Sample initial points
        bounds = torch.tensor(
            [func.lower_bounds, func.upper_bounds], dtype=torch.float32
        )
        X = (torch.rand(num_initial_points, func.dimension)) * (
            bounds[1] - bounds[0]
        ) + bounds[0]
        Y = torch.tensor(
            [func(x.tolist())[0][0] for x in X], dtype=torch.float32
        ).unsqueeze(-1)

        X = normalize(X, bounds=bounds)
        Y = standardize(Y)

        # Fit GP model
        model = SingleTaskGP(X, Y)
        mll = ExactMarginalLogLikelihood(model.likelihood, model)
        fit_gpytorch_mll(mll)
        return model

    @staticmethod
    def update_model(
        model: SingleTaskGP,
        func: StaticProblem,
        new_x: List[List[float]],
        mock: bool = False,
    ) -> SingleTaskGP:
        """
        Update the given GP model with new observations from the function at the specified input points.
        Args:
            model: The existing SingleTaskGP model to be updated.
            func: The static problem function to evaluate at the new input points.
            new_x: A list of new input points (as lists of floats) where the function will be evaluated and the model will be updated.
            mock: If True, only condition the existing model on the new observations without refitting. If False, refit a new model with the combined data.
        Returns:
            An updated SingleTaskGP model that incorporates the new observations from the specified input points.
        """
        new_x_tensor = torch.tensor(new_x, dtype=torch.float32).unsqueeze(0)
        new_y_tensor = torch.tensor(
            [func(x.tolist())[0][0] for x in new_x_tensor], dtype=torch.float32
        ).unsqueeze(-1)
        # Normalize new data point
        bounds = torch.tensor(
            [func.lower_bounds, func.upper_bounds], dtype=torch.float32
        )
        new_x_normalized = normalize(new_x_tensor, bounds=bounds)

        if mock:
            # If mock is True, we only do conditioning
            updated_model = model.condition_on_observations(
                new_x_normalized, new_y_tensor
            )
        else:
            # Combine existing training data with new data
            combined_x = torch.cat([model.train_inputs[0], new_x_normalized], dim=0)
            combined_y = torch.cat([model.train_targets, new_y_tensor], dim=0)

            # Fit a new GP model with the combined data
            updated_model = SingleTaskGP(combined_x, combined_y)
            mll = ExactMarginalLogLikelihood(updated_model.likelihood, updated_model)
            fit_gpytorch_mll(mll)
        return updated_model

    def create_candidate_set(self, model: SingleTaskGP) -> torch.Tensor:
        """
        Create a candidate set for the acquisition function optimization by combining uniformly sampled points and points selected based on the Upper Confidence Bound (UCB) acquisition function.
        Args:
            model: The SingleTaskGP model used to evaluate the UCB acquisition function for selecting candidate points.
        Returns:
            A tensor containing the combined candidate points for acquisition function optimization.
        """
        bounds = torch.tensor(
            [self.func.lower_bounds, self.func.upper_bounds], dtype=torch.float32
        )
        # Sample candidate points uniformly within the bounds
        uniform_candidates = (
            torch.rand(self.num_uniform_candidates, self.func.dimension)
        ) * (bounds[1] - bounds[0]) + bounds[0]
        # Add candidate points based on UCB acquisition function
        ucb_acq_func = UpperConfidenceBound(model, beta=self.ucb_params["beta"])
        ucb_candidates = optimize_acqf(
            acq_function=ucb_acq_func,
            bounds=bounds,
            q=self.num_ucb_candidates,
            num_restarts=self.ucb_params["num_restarts"],
            raw_samples=self.ucb_params["raw_samples"],
        )[0]
        candidate_set = torch.cat([uniform_candidates, ucb_candidates], dim=0)
        return candidate_set

    @staticmethod
    def get_score_func(
        model: SingleTaskGP,
        candidate_set: torch.Tensor,
        acq_func: AcqFunc,
        acq_func_params: Optional[dict] = None,
    ) -> Callable[[List[float]], float]:
        """
        Create a score function based on the specified acquisition function that can be used to evaluate the acquisition value at any given input point.
        Args:
            model: The SingleTaskGP model used to evaluate the acquisition function.
            candidate_set: A tensor containing the candidate points for acquisition function evaluation.
            acq_func: The type of acquisition function to use for scoring (e.g., MAX_VALUE_ENTROPY).
            acq_func_params: Optional parameters specific to the chosen acquisition function.
        Returns:
            A callable function that takes a list of floats as input and returns a float representing the acquisition value at that input point.
        """
        acq_func_params = acq_func_params or {}
        if acq_func == AcqFunc.MAX_VALUE_ENTROPY:

            acquisition_func = qMaxValueEntropy(
                model=model,
                candidate_set=candidate_set,
                **acq_func_params,
            )
        else:
            raise NotImplementedError(
                f"Acquisition function {acq_func} not implemented"
            )

        def score_func(x: List[float]) -> float:
            X = torch.tensor(x, dtype=torch.float32)
            return acquisition_func(X.unsqueeze(0)).item()

        return score_func

    def get_cond_score_func(
        self, train_x: List[List[float]]
    ) -> Callable[[List[float]], float]:
        """
        Get a score function that is conditioned on new observations at the specified input points.
        Args:
                train_x: A list of new input points (as lists of floats) where the function will be evaluated and the model will be updated before creating the score function.
        Returns:
            A callable function that takes a list of floats as input and returns a float representing the acquisition value at that input point, based on the model updated with the new observations.
        """
        # update model with train_x
        tmp_model = ModelImprovementEvaluator.update_model(
            self.model, self.func, train_x
        )
        # TODO: check pass by value etc to ensure self.model is not updated
        candidate_set = self.create_candidate_set(tmp_model)
        score_func = self.get_score_func(
            tmp_model, candidate_set, self.acq_func, self.acq_func_params
        )
        return score_func
