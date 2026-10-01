from abc import ABC, abstractmethod
from os import stat
from pathlib import Path
from mating_kernel.exp.batcher import Batcher
from mating_kernel.problems.mk_suite import MKSuiteConfig
from ttex.config import Config, ConfigurableObject
from mating_kernel.exp.utils.batch import Batch
from mating_kernel.exp.utils.parse_args import experiment_parser
from mating_kernel.exp.factory import parse_experiment_config
import argparse


class ExperimentConfig(Config):
    def __init__(
        self,
        suite_config: MKSuiteConfig,
        settings: dict[str, list],
        reps: int = 1,
        num_batches: int | None = None,
        seed: int | None = None,
        out_dir: str | Path = Path("."),
    ):
        super().__init__()
        self.suite_config = suite_config
        self.reps = reps
        self.settings = settings
        self.num_batches = num_batches
        self.seed = seed
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)


class Experiment(ABC):

    @staticmethod
    @abstractmethod
    def _run_batch(batch: Batch, **kwargs) -> list[str]: ...

    @classmethod
    @abstractmethod
    def parser(cls) -> argparse.ArgumentParser: ...

    @classmethod
    def parse_args(
        cls, argv: list[str] | None = None
    ) -> tuple[ExperimentConfig, list[int]]:
        """
        Create an ExperimentConfig from command line arguments.
        """
        exp_args, unknown_args = experiment_parser().parse_known_args(argv)

        setting_args = cls.parser().parse_args(unknown_args)
        # Create a dictionary of settings from the parsed arguments
        settings_dict = vars(setting_args)

        exp_config = parse_experiment_config(exp_args, settings_dict)
        return exp_config, exp_args.nth

    @classmethod
    def run_from_args(cls, argv: list[str] | None = None, **kwargs):
        """
        Run the experiment from command line arguments.
        """
        config, nth = cls.parse_args(argv)
        return cls.run(config, nth=nth, **kwargs)

    @classmethod
    def run(cls, config: ExperimentConfig, nth: list[int] | None = None, **kwargs):
        """
        Run all batches of the experiment.
        """
        batches = Batcher.create_batches(
            suite_config=config.suite_config,
            settings=config.settings,
            reps=config.reps,
            num_batches=config.num_batches,
            seed=config.seed,
            exp_dir=config.out_dir,
        )
        assert all(
            nth < len(batches) for nth in (nth or [])
        ), "nth indices must be less than the number of batches."
        assert all(nth >= 0 for nth in (nth or [])), "nth indices must be non-negative."
        if nth is not None:
            batches = [b for n in nth for b in [batches[n]]]

        results = []
        for batch in batches:
            result = cls._run_batch(batch, **kwargs)
            results.append(result)
        return results
