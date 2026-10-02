import math
from collections import defaultdict
from itertools import product, zip_longest
from pathlib import Path

import numpy as np

from mating_kernel.exp.utils.batch import Batch
from mating_kernel.problems.mk_suite import MKSuite, MKSuiteConfig


class Batcher:

    @staticmethod
    def seed_batches(
        batches: list[Batch],
        reps: int,
        seed: int | None = None,
    ) -> list[Batch]:
        rng = np.random.default_rng(seed)
        seeds = rng.integers(low=0, high=2**32 - 1, size=reps)
        seeded_batches_dict = defaultdict(list)
        for bid, batch in enumerate(batches):
            existing_seeds = batch.get_run_seeds()
            missing_seed_idx = [
                i for i, s in enumerate(seeds) if s not in existing_seeds
            ]
            for sid in missing_seed_idx:
                seeded_batch = batch.model_copy(
                    update={"seed": int(seeds[sid]), "rep": sid}
                )
                seeded_batches_dict[bid].append(seeded_batch)
        # Flatten the list of seeded batches
        # But intersperse the batch ids
        seeded_batches = [
            batch
            for bid in zip_longest(*seeded_batches_dict.values())
            for batch in bid
            if batch is not None
        ]
        return seeded_batches

    @staticmethod
    def split_batches(
        batches: list[Batch], num_batches: int | None = None
    ) -> list[list[Batch]]:

        if num_batches is None:
            num_batches = len(
                batches
            )  # Default to one batch per batch if not specified
        assert num_batches > 0, "Number of batches must be greater than 0."
        # Split the batches into the specified number of batches
        if num_batches >= len(batches):
            # If the number of requested batches is greater than or equal to the number of batches,
            # return each batch in its own list.
            batched_batches = [[batch] for batch in batches]
        elif num_batches > 1:
            batch_size = math.ceil(len(batches) / num_batches)
            batched_batches = [
                batches[i : i + batch_size] for i in range(0, len(batches), batch_size)
            ]
        else:
            batched_batches = [
                batches
            ]  # Wrap in a list to maintain consistency in return type
        return batched_batches

    @staticmethod
    def create_combinations(
        suite_config: MKSuiteConfig,
        settings: dict[str, list],
        exp_dir: Path | str = Path("."),
    ) -> list[Batch]:
        suite = MKSuite(suite_config)
        problem_ids = list(suite.problem_id_map.keys())
        param_names = list(settings.keys())
        param_values = list(settings.values())
        param_combinations = list(product(*param_values))

        # Collect all combinations of problem IDs and parameter combinations
        # This will create a list of dictionaries, each containing a unique combination of problem ID and parameter settings.
        batches = []
        for problem_id in problem_ids:
            for sid, param_combination in enumerate(param_combinations):
                batch_settings = dict(zip(param_names, param_combination))
                pdata = suite.problem_id_map[problem_id]
                batch = Batch(
                    **batch_settings,
                    pid=problem_id,
                    sid=sid,
                    parent_dir=exp_dir,
                    pinfo=pdata["info"],
                    problem=pdata["problem"],
                )
                batches.append(batch)

        return batches

    @staticmethod
    def create_batches(
        suite_config: MKSuiteConfig,
        settings: dict[str, list],
        reps: int = 1,
        num_batches: int | None = None,  # separate batches if not set
        seed: int | None = None,
        exp_dir: Path | str = Path("."),
        filter_bids: (
            list[int] | None
        ) = None,  # indices of batches to run, if None, run all
    ) -> list[list[Batch]]:
        batches = Batcher.create_combinations(suite_config, settings, exp_dir=exp_dir)
        seeded_batches = Batcher.seed_batches(batches, reps=reps, seed=seed)
        batched_batches = Batcher.split_batches(seeded_batches, num_batches=num_batches)
        if filter_bids is not None:
            if any(bid < 0 or bid >= len(batched_batches) for bid in filter_bids):
                raise ValueError(
                    f"filter_bids contains indices that are out of range. "
                    f"Max index is {len(batched_batches) - 1}."
                )
            # Filter the batched_batches to only include the specified batch ids
            batched_batches = [batched_batches[bid] for bid in filter_bids]

        return batched_batches
