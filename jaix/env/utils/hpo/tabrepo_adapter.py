# from tabrepo.constants.model_constants import MODEL_TYPE_DICT
from tabrepo.repository.evaluation_repository import EvaluationRepository
from enum import Enum
import re
from typing import Optional
import pandas as pd
import numpy as np

# model_types = list(MODEL_TYPE_DICT.values())
# model_tuples = [(mt, mt) for mt in model_types]
# ModelType = Enum("ModelType", model_types)
TaskType = Enum(
    "TaskType", [("C1", "binary"), ("R", "regression"), ("CM", "multiclass")]
)


class TabrepoAdapter:
    def __init__(
        self,
        repo: EvaluationRepository,
        task_type: TaskType,
        inst: int,
    ):
        self.repo = repo
        datasets = self.repo.datasets(union=True, problem_type=task_type.value)
        if inst < len(datasets):
            self.dataset = datasets[inst]
        else:
            raise ValueError(
                f"Tried getting dataset {inst} of {len(datasets)} available"
            )

        # Get all handmade configs that are available for all of the datasets
        regex = r"_c\d+_"
        self.configs = [
            config_name
            for config_name in self.repo.configs(datasets=[self.dataset], union=False)
            if re.search(regex, config_name) is not None
        ]
        self._set_metadata()

    def _set_metadata(self):
        self.metadata = self.repo.dataset_metadata(self.dataset)
        self.metadata.update(self.repo.dataset_info(self.dataset))
        self.metrics = self.repo.metrics(datasets=[self.dataset], configs=self.configs)
        self.metadata["max_error_val"] = max(self.metrics["metric_error_val"])
        self.metadata["min_error_val"] = min(self.metrics["metric_error_val"])

    def evaluate(self, config_id: int, seed: Optional[int] = None):
        # TODO: seed
        idx = pd.IndexSlice
        results = self.metrics.loc[
            idx[:, :, self.configs[config_id]], ["metric_error_val", "time_train_s"]
        ]
        num_folds = len(results.index)
        w = np.random.uniform(size=num_folds)
        w /= sum(w)  # Normalise to act as weights
        metric_error_val = sum(w * results["metric_error_val"])
        time_train_s = sum(w * results["time_train_s"])
        return metric_error_val, time_train_s

    def __str__(self):
        return f"{self.dataset}"
