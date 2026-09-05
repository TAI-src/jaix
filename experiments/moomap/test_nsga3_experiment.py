import numpy as np
import pytest
from jaix.env.utils.archive.mo_archive import MOArchiveConfig
from jaix.env.utils.problem.cobi_problem import CobiProblem
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)

from nsga3_experiment import NSGA3Experiment, NSGA3ExperimentConfig


def test_create_mo_archive_config():
    problem = REProblem(REProblemConfig(), inst=2)
    config = NSGA3ExperimentConfig.create_mo_archive_config(problem)
    assert config.max_size is not None
    assert isinstance(config.num_refpoints, int)
    assert isinstance(config, MOArchiveConfig)


def test_re_problem_list():
    problems = NSGA3ExperimentConfig.re_problem_list()
    assert isinstance(problems, list)
    assert len(problems) > 0
    assert all(isinstance(p, REProblem) for p in problems)


def test_cobi_problem_list():
    problems = NSGA3ExperimentConfig.cobi_problem_list()
    assert isinstance(problems, list)
    assert len(problems) > 0
    assert all(isinstance(p, CobiProblem) for p in problems)


@pytest.mark.parametrize("mode", ["cobi", "reproblem", ""])
def test_generate_problem_list(mode):
    problems = NSGA3ExperimentConfig.generate_problem_list(mode)
    assert isinstance(problems, list)
    assert len(problems) > 0
    if mode == "cobi":
        assert all(isinstance(p, CobiProblem) for p in problems)
    elif mode == "reproblem":
        assert all(isinstance(p, REProblem) for p in problems)
    else:
        assert all(isinstance(p, (CobiProblem, REProblem)) for p in problems)


@pytest.mark.parametrize("kwargs", [{"max_size": 10}, None, {"num_refpoints": 50}])
def test_config_init(kwargs):
    config = NSGA3ExperimentConfig(
        num_independent_runs=1,
        num_generations=10,
        mo_archive_kwargs=kwargs,
    )
    assert config.num_prefill_samples == -1
    problem = REProblem(REProblemConfig(), inst=2)
    config.update_defaults(problem)
    assert isinstance(config.mo_archive_config, MOArchiveConfig)
    assert config.mo_archive_config.max_size is not None
    if kwargs is not None and "max_size" in kwargs:
        assert config.mo_archive_config.max_size == kwargs["max_size"]
    if kwargs is not None and "num_refpoints" in kwargs:
        assert config.mo_archive_config.num_refpoints == kwargs["num_refpoints"]
    assert config.num_prefill_samples == config.mo_archive_config.num_refpoints
    assert config.num_offspring == config.mo_archive_config.num_refpoints


def test_prefill_archive():
    config = NSGA3ExperimentConfig(
        num_independent_runs=1,
        num_generations=1,
        mo_archive_kwargs={"max_size": 10},
    )
    problem = REProblem(REProblemConfig(), inst=2)
    config.update_defaults(problem)
    config.rng = np.random.default_rng(42)
    assert config.num_prefill_samples == config.mo_archive_config.num_refpoints
    archive, entries = NSGA3Experiment.prefill_archive(problem, config)
    assert len(entries) == config.num_prefill_samples
    assert archive.size == config.mo_archive_kwargs["max_size"]

    # Test seeding
    config.rng = np.random.default_rng(42)
    archive2, entries2 = NSGA3Experiment.prefill_archive(problem, config)
    assert len(entries2) == config.num_prefill_samples
    for e1, e2 in zip(entries, entries2):
        assert np.allclose(e1.x, e2.x)
        assert np.allclose(e1.y, e2.y)
    config.rng = np.random.default_rng(43)
    archive3, entries3 = NSGA3Experiment.prefill_archive(problem, config)
    assert len(entries3) == config.num_prefill_samples
    for e1, e3 in zip(entries, entries3):
        assert not np.allclose(e1.x, e3.x) or not np.allclose(e1.y, e3.y)


def test_get_unbounded_archive():
    config = NSGA3ExperimentConfig(
        num_independent_runs=1,
        num_generations=1,
        mo_archive_kwargs={"max_size": 10, "num_refpoints": 5},
    )
    problem = REProblem(REProblemConfig(), inst=2)
    config.update_defaults(problem)
    archive = NSGA3Experiment.get_unbounded_archive(problem, config)
    assert archive.max_size is None
    assert archive.num_refpoints == config.mo_archive_kwargs["num_refpoints"]


def test_entry_info():
    problem = REProblem(REProblemConfig(), inst=2)
    config = NSGA3ExperimentConfig(
        num_independent_runs=1,
        num_generations=1,
        num_prefill_samples=20,
        mo_archive_kwargs={"max_size": 10, "num_refpoints": 5},
    )
    config.update_defaults(problem)
    archive, _entries = NSGA3Experiment.prefill_archive(problem, config)
    for entry in archive.get_all():
        info = NSGA3Experiment.entry_info(archive, entry)
        assert isinstance(info, dict)
        assert "y" in info
        assert "rank" in info
        assert "niche" in info
        assert isinstance(info["rank"], int)
        assert isinstance(info["niche"], int)


def test_select_parents():
    problem = REProblem(REProblemConfig(), inst=2)
    config = NSGA3ExperimentConfig(
        num_independent_runs=1,
        num_generations=1,
        mo_archive_kwargs={"max_size": 10, "num_refpoints": 5},
    )
    config.update_defaults(problem)
    archive, _entries = NSGA3Experiment.prefill_archive(problem, config)
    parents = NSGA3Experiment.select_parents(config, archive)
    assert len(parents) == config.num_parents
    for p_info in parents:
        assert isinstance(p_info, dict)
        assert "rank" in p_info
        assert "niche" in p_info
        assert isinstance(p_info["rank"], int)
        assert isinstance(p_info["niche"], int)
        assert "x" in p_info
        assert isinstance(p_info["x"], (list, tuple, np.ndarray))

    # Test seeding
    config.rng = np.random.default_rng(42)
    parents = NSGA3Experiment.select_parents(config, archive)
    config.rng = np.random.default_rng(42)
    parents2 = NSGA3Experiment.select_parents(config, archive)
    for p1, p2 in zip(parents, parents2):
        assert p1["rank"] == p2["rank"]
        assert p1["niche"] == p2["niche"]
        assert np.allclose(p1["x"], p2["x"])
    config.rng = np.random.default_rng(44)
    parents3 = NSGA3Experiment.select_parents(config, archive)
    for p1, p3 in zip(parents, parents3):
        assert (
            p1["rank"] != p3["rank"]
            or p1["niche"] != p3["niche"]
            or not np.allclose(p1["x"], p3["x"])
        )


def test_create_offspring():
    problem = REProblem(REProblemConfig(), inst=2)
    config = NSGA3ExperimentConfig(
        num_independent_runs=1,
        num_generations=1,
        mo_archive_kwargs={"max_size": 10, "num_refpoints": 5},
    )
    config.update_defaults(problem)
    parents = [
        np.random.uniform(low=problem.lower_bounds, high=problem.upper_bounds)
        for _ in range(config.num_parents)
    ]
    config.rng = np.random.default_rng(42)
    offspring = NSGA3Experiment.create_offspring(parents, config)
    assert len(offspring) == len(parents[0])

    # check that the seeding worked
    config.rng = np.random.default_rng(42)
    offspring2 = NSGA3Experiment.create_offspring(parents, config)
    assert np.allclose(offspring, offspring2)

    # Check a different seed produces different offspring
    config.rng = np.random.default_rng(43)
    offspring3 = NSGA3Experiment.create_offspring(parents, config)
    assert not np.allclose(offspring, offspring3)


def test_create_families():
    problem = REProblem(REProblemConfig(), inst=2)
    config = NSGA3ExperimentConfig(
        num_independent_runs=1,
        num_generations=1,
        num_prefill_samples=20,
        mo_archive_kwargs={"max_size": 5, "num_refpoints": 5},
    )
    config.update_defaults(problem)
    archive, _entries = NSGA3Experiment.prefill_archive(problem, config)
    families = NSGA3Experiment.create_families(config, archive)
    assert len(families) == config.num_offspring
    for family in families:
        assert "parent_0" in family
        assert "offspring" in family
        assert isinstance(family["parent_0"]["rank"], int)
        assert isinstance(family["offspring"], np.ndarray)


def test_try_add_offspring():
    problem = REProblem(REProblemConfig(), inst=2)
    config = NSGA3ExperimentConfig(
        num_independent_runs=1,
        num_generations=1,
        num_prefill_samples=20,
        mo_archive_kwargs={"max_size": 5, "num_refpoints": 5},
    )
    config.update_defaults(problem)
    archive, _entries = NSGA3Experiment.prefill_archive(problem, config)
    offspring = []
    for _ in range(config.num_offspring * 10):
        x = np.random.uniform(low=problem.lower_bounds, high=problem.upper_bounds)
        offspring.append(x)
    offspring_info, _entries = NSGA3Experiment.try_add_offspring(
        problem, archive, offspring
    )
    added = [info["added"] for info in offspring_info]
    assert sum(added) <= archive.max_size
    assert archive.size <= archive.max_size


@pytest.mark.parametrize("mode", ["cobi", "reproblem"])
def test_run_single(tmp_path, mode):
    config = NSGA3ExperimentConfig(
        num_independent_runs=1,
        num_generations=1,
        mode=mode,
        seed=42,
    )
    xp_path = tmp_path / "x1"
    files = NSGA3Experiment.run_single(config, xp_path, problem_idx=None)
    assert len(files) == 2 * len(NSGA3ExperimentConfig.generate_problem_list(mode))
    assert "results" in files[0] and files[0].endswith(".csv")

    # Looking at file 10 since it is the first reproblem (instance 5) with 3 objectives, which means the propulation size is not 100
    with open(files[10], "r") as f:
        import csv

        reader = csv.DictReader(f)
        data = list(reader)
        last_data = data[-1]
        num_data = len(data)
        if mode == "cobi":
            assert num_data == 100  # Cobi problems have dimension 2
        elif mode == "reproblem":
            assert (
                num_data == 91
            )  # reproblem 5 has 3 objectives, so the population size is 91
    assert last_data["seed"] == str(config.seed)

    assert "config" in files[1] and files[1].endswith(".json")
    with open(files[11], "r") as f:
        import json

        config_data = json.load(f)

    assert config_data["NSGA3ExperimentConfig"]["seed"] == config.seed

    # check that seeding worked
    xp_path = tmp_path / "x2"
    # Only run the same problem (instance 5) to check that the results are the same
    files2 = NSGA3Experiment.run_single(config, xp_path, problem_idx=[5])
    with open(files2[0], "r") as f:
        import csv

        reader = csv.DictReader(f)
        data2 = list(reader)
        last_data2 = data2[-1]
    assert last_data2["seed"] == str(config.seed)

    comparison_columns = ["parent_0_x", "parent_1_x", "offspring_x"]
    for col in comparison_columns:
        assert last_data[col] == last_data2[col], f"Mismatch in column {col}"


def test_run_experiment(tmp_path):
    config = NSGA3ExperimentConfig(
        num_independent_runs=2,
        num_generations=2,
        num_prefill_samples=20,
        mode="reproblem",
        seed=42,
    )
    problems_idx = [0, 2]
    xp_path = tmp_path / "experiment"
    config_dict = NSGA3Experiment.run(config, xp_path, problem_idx=problems_idx)
    file = xp_path / f"x_{config_dict['exp_id']}" / "config.json"
    assert file.exists()
    out_files = config_dict["out_files"]
    assert len(out_files) == config.num_independent_runs * 2 * len(problems_idx)

    config_run1 = out_files[1]
    with open(config_run1, "r") as f:
        import json

        config_data = json.load(f)
    seed = config_data["NSGA3ExperimentConfig"]["seed"]
    assert (
        seed is not None and seed != config.seed
    ), f"Seed in config {seed} should not be the same as experiment seed {config.seed}"
