from pathlib import Path
from mating_kernel.exp.utils.find_files import find_data_files


def test_find_data_files():
    # Test that the find_data_files function returns the correct files
    test_folder = Path(__file__).parent.parent.parent / "data" / "nsga3x_results"
    result_files = find_data_files(
        test_folder,
        file_type_pattern="results_*.csv",
        problem_names=["cobi_lin", "cobi_cvex", "test"],
    )
    config_files = find_data_files(test_folder, file_type_pattern="config_*.json")
    assert isinstance(result_files, dict)
    assert isinstance(config_files, dict)
    for files in result_files.values():
        assert all(f.suffix == ".csv" for f in files)
    for files in config_files.values():
        assert all(f.suffix == ".json" for f in files)
    assert set(result_files.keys()) == {"cobi_lin", "cobi_cvex"}
    assert len(result_files["cobi_lin"]) == 2
    assert len(result_files["cobi_cvex"]) == 1
    assert set(config_files.keys()) == {"unknown"}
    result_files2 = find_data_files(
        test_folder, file_type_pattern="results_*.csv", problem_names=["test"]
    )
    assert set(result_files2.keys()) == {"unknown"}
    config_file2 = find_data_files(
        test_folder,
        file_type_pattern="config_*.json",
        problem_names=["cobi_lin", "cobi_cvex", "test"],
    )
    # Check that we have the same order
    for key, cfiles in config_file2.items():
        rfiles = result_files.get(key, [])
        assert len(cfiles) == len(rfiles)
        for cfile, rfile in zip(cfiles, rfiles):
            assert cfile.parent == rfile.parent
            assert cfile.stem.split("_")[-1] == rfile.stem.split("_")[-1]
