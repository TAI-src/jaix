from jaix.env.utils.hpo import ModelType, ConfigTranslator


def test_init():
    model_type = ModelType.XGB
    trl = ConfigTranslator(model_type)
    # https://github.com/autogluon/tabrepo/blob/main/tabrepo/models/xgboost/generate.py
    return trl
