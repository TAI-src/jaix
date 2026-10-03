import argparse
from abc import ABC
import logging

from mating_kernel.exp.utils.batch import Batch
from mating_kernel.exp.utils.batcher import Batcher
from mating_kernel.exp.utils.factory import (
    ExperimentConfig,
    parse_experiment_config,
    ExperimentMode,
)
from mating_kernel.exp.utils.parse_args import experiment_parser

logger = logging.getLogger(__name__)


class Experiment(ABC):
    # TODO: Implement subclass checking to ensure that subclasses implement the required methods.

    @staticmethod
    def _run_batch(batch: Batch, **kwargs) -> list[str]:
        """
        Run a single batch of the experiment.
        """
        raise NotImplementedError(
            "Experiment cls must implement the _run_batch() method to run a single batch."
        )

    @staticmethod
    def parser() -> argparse.ArgumentParser:
        """
        Return an argparse.ArgumentParser for the experiment.
        """
        raise NotImplementedError(
            "Experiment cls must implement the parser() method to return an argparse.ArgumentParser."
        )

    @staticmethod
    def _check_batch_out(batch: Batch, **kwargs) -> bool:
        """
        Check if the batch output is valid.
        """
        raise NotImplementedError(
            "Experiment cls must implement the _check_batch_out() method to check if the batch output is valid."
        )

    @staticmethod
    def _post_process_batch(batch: Batch, **kwargs) -> list[str]:
        """
        Post-process the batch output.
        """
        raise NotImplementedError(
            "Experiment cls must implement the _post_process_batch() method to post-process the batch output."
        )

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
    ) -> list[list[str] | bool]:
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
            skip_existing=(config.mode == ExperimentMode.RUN),
        )
        flattened_batches = [batch for bgroup in batches for batch in bgroup]
        logger.info(
            f"Running Experiment {cls.__name__} mode {config.mode} for {len(flattened_batches)} batches in total."
        )
        results: list[list[str] | bool] = []
        for batch in flattened_batches:
            result: list[str] | bool
            logger.debug(f"Running batch {batch.name}...")
            if config.mode == ExperimentMode.CHECK:
                result = cls._check_batch_out(batch, **kwargs)
            elif config.mode == ExperimentMode.PP:
                result = cls._post_process_batch(batch, **kwargs)
            elif config.mode == ExperimentMode.RUN:
                result = cls._run_batch(batch, **kwargs)
            else:
                raise ValueError(f"Unknown mode: {config.mode}")
            logger.debug(f"Batch {batch.name} result: {result}")
            results.append(result)

        if config.mode == ExperimentMode.CHECK:
            if all(results):
                logger.info("All batches passed the check.")
            else:
                for batch, result in zip(flattened_batches, results):
                    if not result:
                        logger.warning(f"Batch {batch.name} failed the check.")
        else:
            logger.info(f"Finished running {len(results)} batches.")
        return results
