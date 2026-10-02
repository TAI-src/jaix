import argparse
from abc import ABC, abstractmethod

from mating_kernel.exp.utils.batch import Batch
from mating_kernel.exp.utils.batcher import Batcher
from mating_kernel.exp.utils.factory import ExperimentConfig, parse_experiment_config
from mating_kernel.exp.utils.parse_args import experiment_parser


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
        # turn all single values into lists
        for key, value in settings_dict.items():
            if not isinstance(value, list):
                settings_dict[key] = [value]

        exp_config = parse_experiment_config(exp_args, settings_dict)
        return exp_config, exp_args.nth

    @classmethod
    def run_from_args(cls, argv: list[str] | None = None, **kwargs):
        """
        Run the experiment from command line arguments.
        """
        config, filter_bids = cls.parse_args(argv)

        return cls.run(config, filter_bids=filter_bids, **kwargs)

    @classmethod
    def run(
        cls, config: ExperimentConfig, filter_bids: list[int] | None = None, **kwargs
    ):
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
            filter_bids=filter_bids,
        )

        results = []
        for bgroup in batches:
            for batch in bgroup:
                result = cls._run_batch(batch, **kwargs)
                results.append(result)
        return results
