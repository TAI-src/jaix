from jaix.env.utils.archive.moomap_archive import MoomapArchive, MoomapArchiveConfig
from experiment import MoomapX, MoomapXConfig
from jaix.env.wrapper.archive_wrapper import (
    ArchiveWrapperConfig,
)
from jaix.env.wrapper.archive_action_wrapper import (
    ArchiveActionWrapperConfig,
)
from jaix.environment_factory import EnvironmentConfig
from jaix.environment_factory import EnvironmentFactory as EF


def get_config():
    config = MoomapXConfig(
        moomap_config=MoomapArchiveConfig(
            np_bin=1,
            coverage_weight=0.5,
            allow_close_elites=True,
            num_refpoints="original",
        ),
        env_budget_multiplier=1,
        num_samples=100,
        num_trials=10,
        mode="reproblem",
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
            assert entry.info["original_action"] is not None
            break
        break  # Only test one environment for speed


def test_run():
    config = get_config()
    MoomapX.run(config)
