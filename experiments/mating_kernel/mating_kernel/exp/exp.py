import argparse
import logging
import pickle
from abc import ABC
from pathlib import Path

from mating_kernel.exp.utils.batch import Batch
from mating_kernel.exp.utils.batcher import Batcher
from mating_kernel.exp.utils.factory import (
    ExperimentConfig,
    ExperimentMode,
    parse_experiment_config,
)
from mating_kernel.exp.utils.parse_args import experiment_parser

logger = logging.getLogger(__name__)


class Experiment(ABC):
    # TODO: Implement subclass checking to ensure that subclasses implement the required methods.

    @staticmethod
    def _run_batches(batches: list[Batch], **kwargs) -> list[list[str]]:
        """
        Run a single batch of the experiment.
        """
        batch_files = [Experiment.file_paths(batch)["batch"] for batch in batches]
        for batch, batch_file in zip(batches, batch_files):
            logger.debug(f"Finishing batch {batch.name} and saving to {batch_file}")
            with open(batch_file, "wb") as f:
                b_copy = batch.model_copy(update={"problem": str(batch.problem)})
                pickle.dump(b_copy, f)
        return [[str(b_file)] for b_file in batch_files]

    @staticmethod
    def file_paths(batch: Batch) -> dict[str, Path]:
        """
        Return a dictionary of file paths for the batch.
        """
        file_paths_dict = {
            "batch": batch.out_dir / f"batch_{batch.name}.pkl",
        }
        return file_paths_dict

    @staticmethod
    def parser() -> argparse.ArgumentParser:
        """
        Return an argparse.ArgumentParser for the experiment.
        """
        raise NotImplementedError(
            "Experiment cls must implement the parser() method to return an argparse.ArgumentParser."
        )

    @classmethod
    def _check_batches(
        cls, batches: list[Batch], log_level: int = 40, **kwargs
    ) -> list[bool]:
        """
        Check if the batch output is valid.
        """
        checked = [False] * len(batches)
        for i, batch in enumerate(batches):
            files_to_check = list(cls.file_paths(batch).values()) + list(
                Experiment.file_paths(batch).values()
            )
            for file_path in files_to_check:
                if not file_path.exists():
                    logger.log(
                        log_level, f"Batch {batch.name} is missing file: {file_path}"
                    )
                    break
            else:  # If all files exist, mark the batch as checked
                checked[i] = True

        return checked

    @staticmethod
    def _post_process_batches(batches: list[Batch], **kwargs) -> list[list[str]]:
        """
        Post-process the batch output.
        """
        raise NotImplementedError(
            "Experiment cls must implement the _post_process_batches() method to post-process the batch output."
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

    @staticmethod
    def create_batches(
        config: ExperimentConfig, filter_bids: list[int] | None = None
    ) -> list[list[Batch]]:
        """
        Create batches for the experiment.
        """
        batches = Batcher.create_batches(
            suite_config=config.suite_config,
            settings=config.settings,
            reps=config.reps,
            num_batches=config.num_batches,
            seed=config.seed,
            exp_dir=config.out_dir,
            filter_bids=filter_bids,
            skip_existing=config.skip_existing and not config.force_recompute,
            group_by=config.group_by,
        )
        return batches

    @classmethod
    def run(
        cls, config: ExperimentConfig, filter_bids: list[int] | None = None, **kwargs
    ) -> list[list[str] | bool]:
        """
        Run all batches of the experiment.
        """
        batch_groups = Experiment.create_batches(config, filter_bids=filter_bids)
        total_batches = sum(len(bgroup) for bgroup in batch_groups)
        logger.info(
            f"Running Experiment {cls.__name__} mode {config.mode} with {total_batches} batches in {len(batch_groups)} groups."
        )
        results: list[list[str] | bool] = []
        for bgroup in batch_groups:
            result: list[list[str]] | list[bool]
            if config.mode == ExperimentMode.CHECK:
                result = Experiment._check_batches(bgroup, **kwargs)
            elif config.mode == ExperimentMode.PP:
                result = cls._post_process_batches(bgroup, **kwargs)
            elif config.mode == ExperimentMode.RUN:
                if not config.force_recompute:
                    result = Experiment._check_batches(
                        bgroup, log_level=logging.DEBUG, **kwargs
                    )
                    # Skip batches that have already been computed and passed the check
                    if any(result):
                        logger.info(
                            f"Skipping {sum(result)} batches that have already been computed and passed the check."
                        )
                    bgrp_to_run = [
                        batch for batch, res in zip(bgroup, result) if not res
                    ]
                else:
                    bgrp_to_run = bgroup
                cls_files = cls._run_batches(bgrp_to_run, **kwargs)
                exp_files = Experiment._run_batches(bgrp_to_run, **kwargs)
                result = [
                    cls_file + exp_file
                    for cls_file, exp_file in zip(cls_files, exp_files)
                ]
            else:
                raise ValueError(f"Unknown mode: {config.mode}")
            logger.debug(f"Batch group result: {result}")
            results.extend(result)

        if config.mode == ExperimentMode.CHECK:
            if all(results):
                logger.info("All batches passed the check.")
            else:
                flattened_batches = [
                    batch for bgroup in batch_groups for batch in bgroup
                ]
                assert all(
                    isinstance(res, bool) for res in results
                ), "Results should be a list of bools in CHECK mode."
                for batch, res in zip(flattened_batches, results):
                    if not res:
                        logger.warning(f"Batch {batch.name} failed the check.")
                # If any batch failed, raise an exception
                raise RuntimeError(
                    "Some batches failed the check. See logs for details."
                )
        else:
            logger.info(f"Finished running {len(results)} batches.")
        return results
