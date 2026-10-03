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
    def _run_batch(batch: Batch, **kwargs) -> list[str]:
        return [f"ran_{batch.name}"]

    @staticmethod
    def _check_batch_out(batch: Batch, **kwargs) -> bool:
        return True

    @staticmethod
    def _post_process_batch(batch: Batch, **kwargs) -> list[str]:
        return [f"pp_{batch.name}"]


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
    assert results == expected_batches


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
    assert results == expected_batches


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
        assert all(isinstance(res, bool) and res for res in results)
    elif mode == ExperimentMode.PP:
        assert all(
            isinstance(res, list) and res[0].startswith("pp_") for res in results
        )
    elif mode == ExperimentMode.RUN:
        assert all(
            isinstance(res, list) and res[0].startswith("ran_") for res in results
        )
