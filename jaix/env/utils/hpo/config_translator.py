from tabrepo.constants.model_constants import MODEL_TYPE_DICT
from tabrepo.repository.evaluation_repository import load_repository
from enum import Enum

model_types = list(MODEL_TYPE_DICT.values())
model_tuples = [(mt, mt) for mt in model_types]
ModelType = Enum("ModelType", model_types)


class ConfigTranslator:
    def __init__(self, model_type: ModelType, context: str = "D244_F3_C1530_30"):
        self.repo = load_repository(context)
        self.datasets = self.repo.datasets(union=True)
        self.configs = [
            config_name
            for config_name, config_type in self.repo.configs_type().items()
            if config_type == model_type.name
        ]
        self.hyperparams = self.repo.configs_hyperparameters(configs=self.configs)
