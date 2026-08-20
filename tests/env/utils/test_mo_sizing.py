from jaix.env.utils.mo_sizing import (
    get_ref_dirs,
    get_num_refpoints,
    cache_dir,
    cache_enabled,
)


def test_get_ref_dirs():
    ref_dirs = get_ref_dirs(3, "energy_small")
    assert ref_dirs.shape[1] == 3, f"Expected 3 objectives, got {ref_dirs.shape[1]}"
    assert (
        ref_dirs.shape[0] == 50
    ), f"Expected 50 reference points, got {ref_dirs.shape[0]}"


def test_get_num_refpoints():
    # Test for original method
    assert get_num_refpoints(3, "original") == 91
    assert get_num_refpoints(5, "original") == 210
    assert get_num_refpoints(8, "original") == 156
    assert get_num_refpoints(10, "original") == 275
    assert get_num_refpoints(15, "original") == 135

    # Test for energy_small method
    assert get_num_refpoints(3, "energy_small") == 50
    assert get_num_refpoints(5, "energy_small") == 100
    assert get_num_refpoints(8, "energy_small") == 200
    assert get_num_refpoints(10, "energy_small") == 300
    assert get_num_refpoints(15, "energy_small") == 300

    # Test for energy_medium method
    assert get_num_refpoints(3, "energy_medium") == 250
    assert get_num_refpoints(5, "energy_medium") == 250
    assert get_num_refpoints(8, "energy_medium") == 500
    assert get_num_refpoints(10, "energy_medium") == 600
    assert get_num_refpoints(15, "energy_medium") == 600

    # Test for energy_large method
    assert get_num_refpoints(3, "energy_large") == 500
    assert get_num_refpoints(5, "energy_large") == 500
    assert get_num_refpoints(8, "energy_large") == 1000
    assert get_num_refpoints(10, "energy_large") == 1000
    assert get_num_refpoints(15, "energy_large") == 1000

    # Test for inbetween values (interpolation)
    assert get_num_refpoints(4, "original") == 150  # Interpolated between 3 and 5
    assert get_num_refpoints(6, "energy_small") == 133  # Interpolated between 5 and 8


def test_env_values():
    assert cache_enabled is True, "Cache should be enabled by default"
    assert (
        cache_dir == "/tmp/jaix-cache"
    ), f"Expected cache_dir to be '/tmp/jaix-cache', got {cache_dir}"
