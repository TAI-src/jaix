from jaix.env.utils.hpo import TaskType, TabrepoAdapter
import pytest
from ...singular.test_hpo_environment import skip_remaining_tests


@pytest.fixture(scope="session")
def repo(skip_remaining_tests):
    from tabrepo.repository.evaluation_repository import load_repository

    repo = load_repository("D244_F3_C1530_30", load_predictions=True, cache=True)
    return repo


def test_tabrepo_calls(repo):
    # Test that tabrepo calls run
    repo.metrics(
        datasets=["Australian"],
        configs=["CatBoost_r22_BAG_L1", "RandomForest_r12_BAG_L1"],
    )
    repo.evaluate_ensemble(
        datasets=["Australian"],
        configs=[
            "CatBoost_r23_BAG_L1",
            "CatBoost_r22_BAG_L1",
            "RandomForest_r12_BAG_L1",
        ],
        backend="native",
    )


def get_stats(repo):
    # Convenience function to get stats for setting config
    configs, datasets = TabrepoAdapter.get_config_names(repo, task_type=TaskType.C1)
    print("num configs", len(configs), "num functions", len(datasets))
    times = []
    for ds in datasets:
        _, metrics = TabrepoAdapter.get_metadata(repo, ds, configs)
        times.extend(metrics["time_train_s"])
        print(ds, "min time", min(metrics["time_train_s"]))
        print(ds, "max time", max(metrics["time_train_s"]))
        print(
            ds, "mean time", sum(metrics["time_train_s"]) / len(metrics["time_train_s"])
        )
    print("min time", min(times))
    print("max time", max(times))
    print("mean time", sum(times) / len(times))


def test_get_dataset_names(repo):
    datasets = TabrepoAdapter.get_dataset_names(repo, task_type=TaskType.C1)
    assert len(datasets) > 0
    assert datasets[0] == "Australian"

    datasets = TabrepoAdapter.get_dataset_names(
        repo, task_type=TaskType.R, datasets_idx=[0]
    )
    assert len(datasets) == 1

    with pytest.raises(IndexError):
        TabrepoAdapter.get_dataset_names(
            repo, task_type=TaskType.CM, datasets_idx=[300]
        )


def test_get_config_names(repo):
    configs, datasets = TabrepoAdapter.get_config_names(repo, TaskType.C1)
    assert len(configs) > 0
    assert len(datasets) > 0

    configs, ds = TabrepoAdapter.get_config_names(
        repo, TaskType.C1, datasets=datasets[0:2]
    )
    assert len(ds) == 2
    assert len(configs) > 0

    configs, datasets = TabrepoAdapter.get_config_names(repo, TaskType.R, datasets=[0])
    assert len(configs) > 0
    assert len(datasets) == 1


@pytest.mark.parametrize("task_type", [TaskType.C1, TaskType.R, TaskType.CM])
def test_data_summary(repo, task_type):
    configs, datasets = TabrepoAdapter.get_config_names(repo, task_type, datasets=[0])
    metadata, metrics = TabrepoAdapter.get_metadata(repo, datasets[0], configs)
    assert metadata["num_folds"] == 3
    assert metadata["problem_type"] == task_type.value
    assert "max_error_val" in metadata
    assert len(metrics) > 0


@pytest.mark.parametrize(
    "task_type,inst", [(TaskType.C1, 0), (TaskType.R, 1), (TaskType.C1, 300)]
)
def test_adapter(repo, task_type, inst):
    datasets = repo.datasets(union=True, problem_type=task_type.value)
    if inst >= len(datasets):
        with pytest.raises(IndexError):
            adapter = TabrepoAdapter(
                repo=repo, task_type=task_type, dataset_idx=inst, fold=0
            )
        pytest.xfail("Instance does not exist")
    # Continue only if instance exists
    adapter = TabrepoAdapter(repo=repo, task_type=task_type, dataset_idx=inst, fold=0)
    assert adapter.metadata["problem_type"] == task_type.value
    assert adapter.dataset == datasets[inst]

    assert len(adapter.configs) > 0
    for config_name in adapter.configs:
        assert "_c1_" in config_name

    rank, time_train_s = adapter.evaluate(3)
    assert rank >= 0

    metrics = repo.metrics(datasets=[datasets[inst]], configs=[adapter.configs[3]])
    assert min(metrics["time_train_s"]) <= time_train_s
    assert max(metrics["time_train_s"]) >= time_train_s

    #
    # Test evaluate ensemble
    rank_e, time_e = adapter.evaluate_ensemble(config_ids=[0, 1, 3])
    assert rank_e >= 0
    assert rank > rank_e
    assert time_e > time_train_s


def test_evaluate_ensemble(repo):
    adapter = TabrepoAdapter(repo=repo, task_type=TaskType.C1, dataset_idx=0, fold=0)
    rank, time = adapter.evaluate_ensemble(config_ids=[])
    assert rank == adapter.max_rank
    assert time == 0

    rank, time = adapter.evaluate_ensemble(config_ids=[0, 1, 2])
    assert rank < adapter.max_rank
    assert time > 0
