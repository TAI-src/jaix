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
