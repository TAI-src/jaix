# from tabrepo.constants.model_constants import MODEL_TYPE_DICT
from tabrepo.repository.evaluation_repository import EvaluationRepository
from enum import Enum
import re
from typing import Optional, List
import pandas as pd
import numpy as np

# model_types = list(MODEL_TYPE_DICT.values())
# model_tuples = [(mt, mt) for mt in model_types]
# ModelType = Enum("ModelType", model_types)
TaskType = Enum(
    "TaskType", [("C1", "binary"), ("R", "regression"), ("CM", "multiclass")]
)


class TabrepoAdapter:
    def data_summary(
        repo: EvaluationRepository, task_type: TaskType, datasets_idx: List[int] = None
    ):
        dataset_names = repo.datasets(union=True, problem_type=task_type.value)
        data_dict = {}
        if datasets_idx is not None:
            # Will throw an error if not available
            datasets = [dataset_names[i] for i in datasets_idx]
        regex = r"_c\d+_"
        for dataset in datasets:
            # Get all handmade configs that are available for all of the datasets
            regex = r"_c\d+_"
            configs = [
                config_name
                for config_name in repo.configs(datasets=[dataset], union=False)
                if re.search(regex, config_name) is not None
            ]
            metadata, metrics = TabrepoAdapter.get_metadata(repo, dataset, configs)
            # get num_folds
            idx = pd.IndexSlice
            results = metrics.loc[
                idx[:, :, configs[0]], ["metric_error_val", "time_train_s"]
            ]
            num_folds = len(results.index)

            data_dict[dataset] = {
                "configs": configs,
                "metadata": metadata,
                "num_folds": num_folds,
            }
        return data_dict

    def get_metadata(repo: EvaluationRepository, dataset: str, configs: List[str]):
        metadata = repo.dataset_metadata(dataset)
        metadata.update(repo.dataset_info(dataset))
        metrics = repo.metrics(datasets=[dataset], configs=configs)
        metadata["max_error_val"] = max(metrics["metric_error_val"])
        metadata["min_error_val"] = min(metrics["metric_error_val"])
        return metadata, metrics

    def __init__(
        self,
        repo: EvaluationRepository,
        task_type: TaskType,
        dataset_idx: int,
        fold: Optional[int] = None,
    ):
        self.repo = repo
        data_dict = TabrepoAdapter.data_summary(self.repo, task_type, [dataset_idx])

        datasets = self.repo.datasets(union=True, problem_type=task_type.value)
        if dataset_idx < len(datasets):
            self.dataset = datasets[dataset_idx]
        else:
            raise ValueError(
                f"Tried getting dataset {dataset_idx} of {len(datasets)} available"
            )
        # Get all handmade configs that are available for all of the datasets
        regex = r"_c\d+_"
        self.configs = [
            config_name
            for config_name in self.repo.configs(datasets=[self.dataset], union=False)
            if re.search(regex, config_name) is not None
        ]
        self._set_metadata()

        # Set fold
        self.fold = fold
        if fold is not None:
            # Check fold is available
            idx = pd.IndexSlice
            results = self.metrics.loc[
                idx[:, :, self.configs[0]], ["metric_error_val", "time_train_s"]
            ]
            num_folds = len(results.index)
            if fold >= num_folds:
                raise ValueError(
                    f"Tried getting fold {fold} of {num_folds} available for dataset {self.dataset}"
                )

    def evaluate(self, config_id: int, seed: Optional[int] = None):
        assert config_id < len(self.configs)
        # TODO: seed
        if self.fold is None:
            idx = pd.IndexSlice
            results = self.metrics.loc[
                idx[:, :, self.configs[config_id]], ["metric_error_val", "time_train_s"]
            ]
            print(results)
            num_folds = len(results.index)
            w = np.random.uniform(size=num_folds)
            w /= sum(w)  # Normalise to act as weights
            metric_error_val = sum(w * results["metric_error_val"])
            time_train_s = sum(w * results["time_train_s"])
        else:
            idx = pd.IndexSlice
            results = self.metrics.loc[
                idx[:, self.fold, self.configs[config_id]],
                ["metric_error_val", "time_train_s"],
            ]
            metric_error_val = results["metric_error_val"].values[0]
            time_train_s = results["time_train_s"].values[0]

        return metric_error_val, time_train_s

    def evaluate_ensemble(self, config_ids: List[int], seed: Optional[int] = None):
        # print(config_ids < len(self.configs))
        configs = [self.configs[config_id] for config_id in config_ids]
        print(configs)
        print(self.dataset)
        df_out, df_ensemble_weights = self.repo.evaluate_ensemble(
            datasets=[self.dataset], configs=configs, rank=False
        )
        print(df_out)
        print(df_ensemble_weights)

    def __str__(self):
        return f"{self.dataset}"
