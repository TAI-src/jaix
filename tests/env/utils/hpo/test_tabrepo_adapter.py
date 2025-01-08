from jaix.env.utils.hpo import TaskType, TabrepoAdapter
from tabrepo.repository.evaluation_repository import load_repository
import pytest


@pytest.fixture(scope="session")
def repo():
    repo = load_repository("D244_F3_C1530_30", load_predictions=True, cache=True)
    return repo


def test_tabrepo_ö(repo):
    print(
        repo.metrics(
            datasets=["Australian"],
            configs=["CatBoost_r22_BAG_L1", "RandomForest_r12_BAG_L1"],
        )
    )
    print(
        repo.evaluate_ensemble(
            datasets=["Australian"],
            configs=["CatBoost_r22_BAG_L1", "RandomForest_r12_BAG_L1"],
            backend="native",
        )
    )


@pytest.mark.parametrize(
    "task_type,inst", [(TaskType.C1, 0), (TaskType.R, 1), (TaskType.C1, 300)]
)
def test_init(repo, task_type, inst):
    datasets = repo.datasets(union=True, problem_type=task_type.value)
    if inst >= len(datasets):
        with pytest.raises(ValueError):
            adapter = TabrepoAdapter(repo=repo, task_type=task_type, dataset_idx=inst)
        pytest.xfail("Instance does not exist")
    # Continue only if instance exists
    adapter = TabrepoAdapter(repo=repo, task_type=task_type, dataset_idx=inst)
    assert adapter.metadata["problem_type"] == task_type.value
    assert adapter.dataset == datasets[inst]

    assert len(adapter.configs) > 0
    hyperparams = repo.configs_hyperparameters(configs=adapter.configs)
    for config_name, params in hyperparams.items():
        assert params["ag_args"]["name_suffix"].startswith("_c")

    metric_error_val, time_train_s = adapter.evaluate(3)
    assert metric_error_val <= adapter.metadata["max_error_val"]
    assert metric_error_val >= adapter.metadata["min_error_val"]

    metrics = repo.metrics(datasets=[datasets[inst]], configs=[adapter.configs[3]])
    assert min(metrics["metric_error_val"]) <= metric_error_val
    assert max(metrics["metric_error_val"]) >= metric_error_val
    assert min(metrics["time_train_s"]) <= time_train_s
    assert max(metrics["time_train_s"]) >= time_train_s

    # TODO: test single fold
    #
    # Test evaluate ensemble
    adapter.evaluate_ensemble(config_ids=[0, 1, 2])
