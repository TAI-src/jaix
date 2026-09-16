from pathlib import Path

from utils_read import (
    find_data_files,
    get_config_dict,
    get_nsga3x_results,
    read_perf_results,
)


def test_get_config_file():
    # Test that the config file is read correctly
    test_config_file = Path(__file__).parent / "data" / "test_config.json"
    config_dict = get_config_dict(test_config_file)
    assert isinstance(config_dict, dict)


def test_find_data_files():
    # Test that the find_data_files function returns the correct files
    test_folder = Path(__file__).parent / "data" / "nsga3x_results"
    result_files = find_data_files(test_folder, file_type_pattern="results_*.csv")
    config_files = find_data_files(test_folder, file_type_pattern="config_*.json")
    assert isinstance(result_files, dict)
    assert isinstance(config_files, dict)
    for files in result_files.values():
        assert all(f.suffix == ".csv" for f in files)
    for files in config_files.values():
        assert all(f.suffix == ".json" for f in files)
    assert set(result_files.keys()) == {0, 1}
    assert len(result_files[0]) == 2
    assert len(result_files[1]) == 1
    for config, problem in zip(config_files[0], result_files[0]):
        assert (
            config.parent.name == problem.parent.name
        ), f"Config file {config} and result file {problem} are not in the same folder"

    result_files2 = find_data_files(
        test_folder, file_type_pattern="results_*.csv", problem_ids=[0]
    )
    assert set(result_files2.keys()) == {0}


def test_get_nsga3x_results():
    # Test that the get_nsga3x_results function returns the correct results
    test_folder = Path(__file__).parent / "data" / "nsga3x_results"
    problem_ids = [0, 1, 5, 7]
    results = get_nsga3x_results(test_folder, problem_ids=problem_ids)
    assert isinstance(results, dict)
    assert set(results.keys()) == set(problem_ids)
    for problem_id, problem_info in results.items():
        assert isinstance(problem_info, dict)
        assert "nadir_point" in problem_info
        run_ids = [rid for rid in problem_info if rid.startswith("r_")]
        if problem_id == 0:
            assert len(run_ids) == 2
        elif problem_id == 1:
            assert len(run_ids) == 1
        else:
            assert len(run_ids) == 0
        for run_id in run_ids:
            run_info = problem_info[run_id]
            assert "result_file" in run_info
            assert "config_file" in run_info
            assert "config" in run_info
            assert isinstance(run_info["config"], dict)


def test_read_perf_results():
    # Test that the read_perf_results function returns the correct results
    data_folder = Path(__file__).parent.parent / "tmp_res"
    read_results = read_perf_results(
        data_folder, algorithm_names=["nsga3"], problem_ids=None
    )
