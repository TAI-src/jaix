from pathlib import Path

from jaix.env.utils.archive.moomap_archive import (
    EliteSelectionStrategy,
    MoomapArchiveConfig,
)
from jaix.env.wrapper.archive_action_wrapper import (
    ArchiveActionWrapperConfig,
)
from jaix.env.wrapper.archive_wrapper import (
    ArchiveWrapperConfig,
)
from jaix.environment_factory import EnvironmentConfig
from jaix.environment_factory import EnvironmentFactory as EF

from experiment import MoomapX, MoomapXConfig, main


def get_config(mode="reproblem") -> MoomapXConfig:
    config = MoomapXConfig(
        moomap_config=MoomapArchiveConfig(
            np_bin=1,
            coverage_weight=0.5,
            allow_close_elites=True,
            num_refpoints="original",
            elite_selection_strategy=EliteSelectionStrategy.FITPROP,
        ),
        num_samples=2,
        num_trials=10,
        mode=mode,
        seed=1337,
    )
    return config


def test_gen_aw_config():
    awc = get_config().generate_archive_wrapper_config()
    assert isinstance(awc, ArchiveWrapperConfig)


def test_gen_aawc_config():
    aawc = get_config().generate_archive_action_wrapper_config()
    assert isinstance(aawc, ArchiveActionWrapperConfig)


def test_gen_env_config():
    env_config = get_config().generate_env_config()
    assert isinstance(env_config, EnvironmentConfig)


def test_prefill():
    config = get_config()
    env_config = config.generate_env_config()
    for env in EF.get_envs(env_config):
        archive = MoomapX.prefill_archive(env, config)
        entries = archive.get_all()
        assert len(entries) > 0
        for entry in entries:
            assert entry.info["env_action"] is not None
            break
        break  # Only test one environment for speed


def test_run(tmp_path):
    config = get_config()
    MoomapX.run(config, out_dir=tmp_path)
    assert (tmp_path / "experiment_info.json").exists()
    # check the number of csv files in the outdir. It should be one for each instance of reproblems
    num_csv_files = len(list(tmp_path.glob("*.csv")))
    assert num_csv_files == 23, f"Expected 23 CSV files, but found {num_csv_files}"


def test_main(tmp_path):
    config_file = "test_config.json"
    out_dir = main(config_file=config_file, out_dir=tmp_path)
    # check that there is a folder with experiment id in the outpath
    assert Path(out_dir).exists() and Path(out_dir).is_dir()
    # check that there is a experiment_info.json file in the outpath
    assert (Path(out_dir) / "experiment_info.json").exists()


def test_cobi(tmp_path):
    try:
        import cobi

        assert cobi is not None
    except ImportError:
        # skipt the test if cobi is not installed
        return
    config = get_config(mode="cobi")
    exp_info = MoomapX.run(config, out_dir=tmp_path)
    assert len(exp_info["envs"]) == 7
