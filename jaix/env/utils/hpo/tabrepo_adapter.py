from tabrepo.constants.model_constants import MODEL_TYPE_DICT
from tabrepo.repository.evaluation_repository import EvaluationRepository
import re
from typing import Optional, List, Union, Tuple
import pandas as pd
import numpy as np
from jaix.env.utils.hpo import TaskType

# model_types = list(MODEL_TYPE_DICT.values())
# model_tuples = [(mt, mt) for mt in model_types]
# ModelType = Enum("ModelType", model_types)


class TabrepoAdapter:
    @staticmethod
    def get_dataset_names(
        repo: EvaluationRepository,
        task_type: TaskType,
        datasets_idx: Optional[List[int]] = None,
    ) -> List[str]:
        dataset_names = repo.datasets(union=True, problem_type=task_type.value)
        if datasets_idx is not None:
            # Will throw an error if not available
            datasets = [dataset_names[i] for i in datasets_idx]
        else:
            datasets = dataset_names
        return datasets

    @staticmethod
    def get_config_names(
        repo: EvaluationRepository,
        task_type: TaskType,
        datasets: Optional[Union[List[int], List[str]]] = None,
    ) -> Tuple[List[str], List[str]]:
        if datasets is None or isinstance(datasets[0], int):
            dataset_names: List[str] = TabrepoAdapter.get_dataset_names(
                repo,
                task_type,
                datasets,  # type: ignore
            )
        else:
            dataset_names: List[str] = datasets  # type: ignore
        # Get one config per type
        regex = r"_c1_"  # For all handmade: r"_c\d+_"
        configs = [
            config_name
            for config_name in repo.configs(datasets=dataset_names, union=False)
            if re.search(regex, config_name) is not None
        ]
        return configs, dataset_names

    @staticmethod
    def get_metadata(repo: EvaluationRepository, dataset: str, configs: List[str]):
        metadata = repo.dataset_metadata(dataset)
        metadata.update(repo.dataset_info(dataset))
        metrics = repo.metrics(datasets=[dataset], configs=configs)
        idx = pd.IndexSlice
        results = metrics.loc[
            idx[:, :, configs[0]], ["metric_error_val", "time_train_s"]
        ]
        metadata["max_error_val"] = max(metrics["metric_error_val"])
        metadata["min_error_val"] = min(metrics["metric_error_val"])
        metadata["mean_training_time"] = sum(results["time_train_s"]) / len(
            results["time_train_s"]
        )

        num_folds = len(results.index)
        metadata["num_folds"] = num_folds
        return metadata, metrics

    def __init__(
        self,
        repo: EvaluationRepository,
        task_type: TaskType,
        dataset_idx: int,
        fold: int,
    ):
        self.repo = repo
        self.configs, datasets = TabrepoAdapter.get_config_names(self.repo, task_type)
        self.dataset = datasets[dataset_idx]
        self.metadata, self.metrics = TabrepoAdapter.get_metadata(
            self.repo, self.dataset, self.configs
        )

        self.max_rank = 1530
        # Set fold
        self.fold = fold
        if fold is not None and fold >= self.metadata["num_folds"]:
            raise ValueError(
                f"Tried getting fold {fold} of {self.metadata['num_folds']} available for dataset {self.dataset}"
            )

    def evaluate(self, config_id: int, seed: Optional[int] = None):
        assert config_id < len(self.configs)
        # TODO: seed
        idx = pd.IndexSlice
        results = self.metrics.loc[
            idx[:, self.fold, self.configs[config_id]],
            ["rank", "time_train_s"],
        ]
        rank = results["rank"].values[0]
        time_train_s = results["time_train_s"].values[0]
        return rank, time_train_s

    def evaluate_ensemble(self, config_ids: List[int], seed: Optional[int] = None):
        if len(config_ids) == 0:
            return self.max_rank, 0
        assert max(config_ids) < len(self.configs)
        configs = [self.configs[config_id] for config_id in config_ids]

        df_out, _ = self.repo.evaluate_ensemble(
            datasets=[self.dataset],
            configs=configs,
            rank=True,
            folds=[self.fold],
            backend="native",  # This makes it faster
        )
        rank = df_out.iloc[0]
        time_train_s = sum(
            [self.evaluate(config_id, seed=seed)[1] for config_id in config_ids]
        )
        return rank, time_train_s

    def __str__(self):
        return f"{self.dataset}"
