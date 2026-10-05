import argparse
from pathlib import Path
from mating_kernel.problems.mk_suite import MKSuiteConfig
import pytest

from mating_kernel.exp.exp import Experiment
from mating_kernel.exp.utils.batch import Batch
from mating_kernel.exp.utils.factory import ExperimentConfig, ExperimentMode


class DummyExperiment(Experiment):
    @staticmethod
    def parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument(
            "--foo",
            type=int,
            nargs="+",
            default=None,
            help="Dummy argument for testing.",
        )
        parser.add_argument(
            "--bar",
            type=str,
            default="default",
            help="Another dummy argument for testing.",
        )
        return parser

    @staticmethod
    def file_paths(batch: Batch) -> dict[str, Path]:
        file_paths_dict = {
            "run": batch.out_dir / f"ran_{batch.name}",
        }
        return file_paths_dict

    @staticmethod
    def _run_batches(batches: list[Batch], **kwargs) -> list[list[str]]:
        res_files = []
        for batch in batches:
            files = list(DummyExperiment.file_paths(batch).values())
            str_files = [str(f) for f in files]
            res_files.append(str_files)
        return res_files

    @staticmethod
    def _post_process_batches(batches: list[Batch], **kwargs) -> list[list[str]]:
        pp_file = batches[0].out_dir / f"pp_{batches[0].name}"
        return [[str(pp_file)] for batch in batches]


def test_parse_args():
    config, filter_bids = DummyExperiment.parse_args(
        [
            "--cobi",
            "--reps",
            "5",
            "--num_batches",
            "2",
            "--seed",
            "123",
            "--out_dir",
            "results",
            "--nth",
            "1",
            "3",
            "--foo",
            "99",
            "100",
            "--group_by",
            "pid",
            "sid",
        ]
    )

    assert config.suite_config.cobi is True
    assert config.reps == 5
    assert config.num_batches == 2
    assert config.seed == 123
    assert config.out_dir == Path("results")

    assert config.settings == {"foo": [99, 100], "bar": ["default"]}
    assert config.group_by == ["pid", "sid"]

    assert filter_bids == [1, 3]


def test_run(tmp_path):
    config, filter_bids = DummyExperiment.parse_args(
        [
            "--cobi",
            "--reps",
            "5",
            "--seed",
            "123",
            "--out_dir",
            str(tmp_path),
            "--foo",
            "99",
            "100",
            "--num_batches",
            "1",
        ]
    )
    assert filter_bids is None

    results = DummyExperiment.run(config, filter_bids=filter_bids)
    from mating_kernel.problems.cobi_configs import names as cobi_names

    expected_batches = []
    for rep in range(config.reps):
        for cobi_name in cobi_names:
            for sid in range(2):  # 2 different settings based on --foo
                expected_batches.append([f"ran_p{cobi_name}_s{sid}_r{rep}"])
    assert len(results) == len(expected_batches)
    res_names = [[Path(r[0]).name] for r in results]
    assert res_names == expected_batches


def test_run_from_args(tmp_path):
    results = DummyExperiment.run_from_args(
        [
            "--cobi",
            "--reps",
            "5",
            "--seed",
            "123",
            "--out_dir",
            str(tmp_path),
            "--foo",
            "99",
            "100",
            "--num_batches",
            "1",
        ]
    )
    from mating_kernel.problems.cobi_configs import names as cobi_names

    expected_batches = []
    for rep in range(5):
        for cobi_name in cobi_names:
            for sid in range(2):  # 2 different settings based on --foo
                expected_batches.append([f"ran_p{cobi_name}_s{sid}_r{rep}"])
    assert len(results) == len(expected_batches)
    res_names = [[Path(r[0]).name] for r in results]
    assert res_names == expected_batches


def test_check_batch_out(tmp_path):
    config = ExperimentConfig(
        suite_config=MKSuiteConfig(
            cobi=True, re=True, constrained=False, num_objectives=[4]
        ),
        settings={"foo": [1, 2], "bar": ["a"]},
        reps=2,
        num_batches=1,
        seed=123,
        out_dir=tmp_path,
        mode=ExperimentMode.CHECK,
    )
    batches = DummyExperiment.create_batches(config)
    for bgroup in batches:
        assert not any(DummyExperiment._check_batches(bgroup))
        # Now create the expected output files for each batch
        run_files = DummyExperiment._run_batches(bgroup)
        flattened_files = [f for sublist in run_files for f in sublist]
        for file_path in flattened_files:
            Path(file_path).parent.mkdir(parents=True, exist_ok=True)
            Path(file_path).touch()
        # We didn't create the batch.pkl files, so the check should still fail
        assert not any(DummyExperiment._check_batches(bgroup))
        batch_pkl_files = Experiment._run_batches(
            bgroup
        )  # This will create the batch.pkl files
        flattened_files = [f for sublist in batch_pkl_files for f in sublist]
        for file_path in flattened_files:
            Path(file_path).parent.mkdir(parents=True, exist_ok=True)
            Path(file_path).touch()

        # Now the check should pass
        assert all(DummyExperiment._check_batches(bgroup))


@pytest.mark.parametrize(
    "mode",
    [ExperimentMode.CHECK, ExperimentMode.PP, ExperimentMode.RUN],
)
def test_run_check_mode(mode, tmp_path):
    config = ExperimentConfig(
        suite_config=MKSuiteConfig(
            cobi=True, re=True, constrained=False, num_objectives=[4]
        ),
        settings={"foo": [1, 2], "bar": ["a"]},
        reps=2,
        num_batches=1,
        seed=123,
        out_dir=tmp_path,
        mode=mode,
    )
    if mode == ExperimentMode.CHECK:
        # We are expecting the check to fail because we haven't created any output files yet
        with pytest.raises(RuntimeError):
            results = DummyExperiment.run(config)
    else:
        results = DummyExperiment.run(config)
        assert len(results) == 2 * 2 * 2  # 2 problems * 2 settings * 2 reps
        if mode == ExperimentMode.PP:
            assert all(isinstance(res, list) and "pp_" in res[0] for res in results)
        elif mode == ExperimentMode.RUN:
            assert all(isinstance(res, list) and "ran_" in res[0] for res in results)
            assert all(isinstance(res, list) and "batch_" in res[1] for res in results)


@pytest.mark.parametrize(
    "skip_existing",
    [True, False],
)
def test_skip_existing(skip_existing, tmp_path):
    config = ExperimentConfig(
        suite_config=MKSuiteConfig(
            cobi=True, re=True, constrained=False, num_objectives=[4]
        ),
        settings={"foo": [1], "bar": ["a"]},
        reps=1,
        num_batches=1,
        seed=123,
        out_dir=tmp_path,
        mode=ExperimentMode.RUN,
        skip_existing=skip_existing,
    )
    # Run the experiment once to create the output files
    results_first_run = DummyExperiment.run(config)
    assert len(results_first_run) == 2  # 2 problems * 1 setting * 1 rep
    # Run the experiment again with skip_existing=True
    # Modify config to create more batches (that have not been run)
    config.reps = 2  # Increase reps to create new batches

    results_second_run = DummyExperiment.run(config)
    # The second run should only process the new batches (the ones that were not run in the first run)
    assert len(results_second_run) == 2  # 2 problems * 1 setting * 1 new rep
    # check that the first batch was skipped and the second batch was run
    for res in results_second_run:
        assert "s0_r1" in res[0]
    # If we run now, everything should be skipped
    results_third_run = DummyExperiment.run(config)
    assert len(results_third_run) == 0


@pytest.mark.parametrize("skip_existing", [True, False])
def test_force_recompute(skip_existing, tmp_path):
    config = ExperimentConfig(
        suite_config=MKSuiteConfig(
            cobi=True, re=True, constrained=False, num_objectives=[4]
        ),
        settings={"foo": [1], "bar": ["a"]},
        reps=1,
        num_batches=1,
        seed=123,
        out_dir=tmp_path,
        mode=ExperimentMode.RUN,
        skip_existing=skip_existing,
        force_recompute=False,
    )
    # Run the experiment once to create the output files
    results_first_run = DummyExperiment.run(config)
    assert len(results_first_run) == 2  # 2 problems * 1 setting * 1 rep
    # Run again, should not have new results
    results_second_run = DummyExperiment.run(config)
    assert len(results_second_run) == 0  # No new results since nothing changed

    # Now set force_recompute to True and run again
    config.force_recompute = True
    results_third_run = DummyExperiment.run(config)
    assert len(results_third_run) == 2  # Should recompute all batches
