import argparse
import copy
import logging
import pickle
from pathlib import Path

from mating_kernel.pymoo.mating.preference_mating import PreferenceMating
from mating_kernel.pymoo.pref_model.pref_model import PreferenceModel
import pandas as pd
from pymoo.algorithms.moo.nsga2 import NSGA2, binary_tournament
from pymoo.core.algorithm import Algorithm
from pymoo.core.result import Result
from pymoo.optimize import minimize

from jaix.env.utils.problem.static_problem import StaticProblem
from mating_kernel.exp.exp import Experiment
from mating_kernel.exp.utils.batch import Batch
from mating_kernel.pymoo.do_recorder import RecordingConfig
from mating_kernel.pymoo.offspring_success_recording_callback import (
    OffspringSuccessRecordingCallback,
)
from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper
from mating_kernel.pymoo.recordable_object import make_recordable
from mating_kernel.exp.utils.get_default import get_default

logger = logging.getLogger(__name__)


class PerfExperiment(Experiment):

    @staticmethod
    def parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument(
            "--alg_name",
            nargs="+",
            type=str,
            default="NSGA2",
            choices=["NSGA2"],
            help="Name of the algorithm to run.",
        )
        parser.add_argument(
            "--pref_model",
            nargs="+",
            type=str,
            default="default",
        )
        parser.add_argument(
            "--n_gen",
            type=int,
            default=1000,
            help="Number of generations to run the algorithm for.",
        )
        parser.add_argument(
            "--num_candidates",
            type=int,
            default=1,
            nargs="+",
            help="Number of candidates for preference selection",
        )
        parser.add_argument(
            "--candidate_pressure",
            type=int,
            default=2,
            nargs="+",
            help="Pressure for candidate selection",
        )
        parser.add_argument(
            "--num_oracle_simulations",
            type=int,
            default=30,
            nargs="+",
            help="Number of times to simulate the offspring for oracle preference selection",
        )

        return parser

    @staticmethod
    def file_paths(batch: Batch) -> dict[str, Path]:
        file_paths_dict = {
            "result": batch.out_dir / f"result_{batch.name}.pkl",
            "record_stats": batch.out_dir / f"record_stats_{batch.name}.csv",
            "archive": batch.out_dir / f"archive_{batch.name}.pkl",
        }
        return file_paths_dict

    @staticmethod
    def _run_batches(batches: list[Batch], **kwargs) -> list[list[str]]:
        batch_files = []
        for b in batches:
            logger.debug(f"Running batch {b.name}")
            files = PerfExperiment._run_batch(b, **kwargs)
            batch_files.append(files)
        return batch_files

    @staticmethod
    def _run_batch(batch: Batch, **kwargs) -> list[str]:
        pref_attrs = [
            "num_oracle_simulations",
        ]

        pref_params = {
            attr: getattr(batch, attr) for attr in pref_attrs if hasattr(batch, attr)
        }
        result, record_stats = PerfExperiment.run_instrumented_pymoo(
            problem=batch.problem,
            algorithm_name=batch.alg_name,  # type: ignore[attr-defined]
            pref_model=batch.pref_model,  # type: ignore[attr-defined]
            pref_params=pref_params,
            num_candidates=batch.num_candidates,  # type: ignore[attr-defined]
            candidate_pressure=batch.candidate_pressure,  # type: ignore[attr-defined]
            n_gen=batch.n_gen,  # type: ignore[attr-defined]
            algorithm_params={},  # type: ignore[attr-defined]
            seed=batch.seed,
        )
        files = PerfExperiment.file_paths(batch)

        # Save the result and record_stats to files in batch.out_dir
        with open(files["result"], "wb") as f:
            result_cpy = copy.deepcopy(result)
            result_cpy.problem = str(batch.problem)
            result_cpy.algorithm = batch.alg_name  # type: ignore[attr-defined]
            pickle.dump(result_cpy, f)
        record_stats.to_csv(files["record_stats"], index=False)
        # Save the archive as pkl as well
        if hasattr(batch.problem, "archive") and hasattr(
            batch.problem.archive, "archived_entries"
        ):
            archive_entries = batch.problem.archive.archived_entries
        else:
            archive_entries = []  # If no archive exists, save an empty list
        with open(files["archive"], "wb") as f:
            pickle.dump(archive_entries, f)

        return [str(f) for f in files.values()]

    @staticmethod
    def _post_process_batches(batches: list[Batch], **kwargs) -> list[list[str]]:
        # Expecting batches to all be from the experiment, i.e. have the same expeirment id. (this means same problem and same settings)
        exp_id = batches[0].experiment_id
        assert all(
            b.experiment_id == exp_id for b in batches
        ), "All batches must have the same experiment id"
        # Merge the record_stats files into one
        record_stats_files = [
            PerfExperiment.file_paths(b)["record_stats"] for b in batches
        ]
        seeds = [b.seed for b in batches]
        merged_record_stats = pd.concat(
            [
                pd.read_csv(f).assign(seed=seed)
                for f, seed in zip(record_stats_files, seeds)
            ],
            ignore_index=True,
        )
        # Save the merged record_stats to a new file
        out_path = Path(batches[0].parent_dir) / f"merged_record_stats_{exp_id}.csv"
        merged_record_stats.to_csv(out_path, index=False)
        return [[str(out_path)]]

    @staticmethod
    def get_alg_class(algorithm_name: str) -> Algorithm:
        if algorithm_name == "NSGA2":
            return NSGA2
        else:
            raise ValueError(f"Unsupported algorithm: {algorithm_name}")

    @staticmethod
    def get_recorded_alg(
        algorithm_name: str,
        pref_model: str | None,
        pref_params: dict | None,
        num_candidates: int,
        candidate_pressure: int,
        algorithm_params: dict | None,
        record_args: dict[str, RecordingConfig] | None = None,
        record_attributes: dict[str, RecordingConfig] | None = None,
    ) -> Algorithm:
        algorithm_class = PerfExperiment.get_alg_class(algorithm_name)
        record_alg_class = make_recordable(algorithm_class)
        algorithm_params = algorithm_params or {}
        pref_params = pref_params or {}
        pref_mod: None | PreferenceModel
        if pref_model in ("default", None):
            pref_mod = None
        elif pref_model == "random":
            from mating_kernel.pymoo.pref_model.random_pref_model import (
                RandomPreferenceModel,
            )

            pref_mod = RandomPreferenceModel(**pref_params)
        elif pref_model == "oracle":
            from mating_kernel.pymoo.pref_model.oracle_pref_model import (
                OraclePreferenceModel,
            )

            # FIXME: This is a bit hacky
            func_comp = pref_params.pop("func_comp", None)
            if func_comp is None:
                func_comp = (
                    binary_tournament  # Default to binary tournament if not provided
                )
            pref_mod = OraclePreferenceModel(func_comp=func_comp, **pref_params)
        else:
            raise ValueError(f"Unsupported preference model: {pref_model}")
        if pref_mod is not None:
            mating_operators = {}
            for attr in ["selection", "crossover", "mutation"]:
                operator = algorithm_params.pop(attr, None)
                if operator is None:
                    operator = get_default(algorithm_class, attr)
                mating_operators[attr] = operator
            mating = PreferenceMating(
                selection=mating_operators["selection"],
                crossover=mating_operators["crossover"],
                mutation=mating_operators["mutation"],
                pref_model=pref_mod,
                num_candidates=num_candidates,
                candidate_pressure=candidate_pressure,
            )
            algorithm_params["mating"] = mating

        algorithm = record_alg_class(
            **algorithm_params,
            record_args=record_args,
            record_attributes=record_attributes,
        )
        return algorithm

    @staticmethod
    def run_instrumented_pymoo(
        problem: StaticProblem,
        algorithm_name: str,
        pref_model: str | None,
        pref_params: dict | None,
        num_candidates: int,
        candidate_pressure: int,
        n_gen: int,
        algorithm_params: dict | None = None,
        seed: int | None = None,
    ) -> tuple[Result, pd.DataFrame]:
        pymoo_problem = PymooProblemWrapper(problem)
        callback = OffspringSuccessRecordingCallback(pymoo_problem)
        algorithm = PerfExperiment.get_recorded_alg(
            algorithm_name=algorithm_name,
            pref_model=pref_model,
            pref_params=pref_params,
            num_candidates=num_candidates,
            candidate_pressure=candidate_pressure,
            algorithm_params=algorithm_params,
            record_args=callback.record_arg_keys,
            record_attributes=callback.record_attribute_keys,
        )

        result = minimize(
            pymoo_problem,
            algorithm,
            seed=seed,
            termination=("n_gen", n_gen),
            callback=callback,
            verbose=logger.isEnabledFor(logging.DEBUG),
        )
        return result, callback.record_stats
