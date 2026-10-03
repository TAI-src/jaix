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
    def _run_batch(batch: Batch, **kwargs) -> list[str]:
        files = list(DummyExperiment.file_paths(batch).values())
        return [str(f) for f in files]

    @staticmethod
    def _post_process_batch(batch: Batch, **kwargs) -> list[str]:
        pp_file = batch.out_dir / f"pp_{batch.name}"
        return [str(pp_file)]


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
        ]
    )

    assert config.suite_config.cobi is True
    assert config.reps == 5
    assert config.num_batches == 2
    assert config.seed == 123
    assert config.out_dir == Path("results")

    assert config.settings == {"foo": [99, 100], "bar": ["default"]}

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
    for batch in batches:
        assert not DummyExperiment._check_batch_out(batch)
        # Now create the expected output files for each batch
        for file_path in DummyExperiment.file_paths(batch).values():
            file_path.parent.mkdir(parents=True, exist_ok=True)
            file_path.touch()  # Create an empty file
        for file_path in Experiment.file_paths(batch).values():
            file_path.parent.mkdir(parents=True, exist_ok=True)
            file_path.touch()  # Create an empty file
        assert DummyExperiment._check_batch_out(batch)


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
    results = DummyExperiment.run(config)
    assert len(results) == 2 * 2 * 2  # 2 problems * 2 settings * 2 reps
    if mode == ExperimentMode.CHECK:
        assert all(isinstance(res, bool) and not res for res in results)
    elif mode == ExperimentMode.PP:
        assert all(isinstance(res, list) and "pp_" in res[0] for res in results)
    elif mode == ExperimentMode.RUN:
        assert all(isinstance(res, list) and "ran_" in res[0] for res in results)
        assert all(isinstance(res, list) and "batch_" in res[1] for res in results)
