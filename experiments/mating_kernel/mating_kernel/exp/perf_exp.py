import argparse
import copy
import logging
import pickle
from pathlib import Path
from typing import ClassVar
from cocoviz import ResultSet

import pandas as pd
from pymoo.algorithms.moo.nsga2 import NSGA2, binary_tournament
from pymoo.core.algorithm import Algorithm
from pymoo.core.result import Result
from pymoo.optimize import minimize

from jaix.env.utils.problem.static_problem import StaticProblem
from mating_kernel.exp.exp import Experiment
from mating_kernel.exp.utils.batch import Batch
from mating_kernel.exp.utils.factory import ExperimentMode, ResultsConfig
from mating_kernel.exp.utils.get_default import get_default
from mating_kernel.pymoo.do_recorder import RecordingConfig
from mating_kernel.pymoo.mating.preference_mating import PreferenceMating
from mating_kernel.pymoo.offspring_success_recording_callback import (
    OffspringSuccessRecordingCallback,
)
from mating_kernel.pymoo.pref_model.pref_model import PreferenceModel
from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper
from mating_kernel.pymoo.recordable_object import make_recordable
from mating_kernel.plots.cv_result import Result as CVResult
from mating_kernel.plots.rt_plot import plot_rt

logger = logging.getLogger(__name__)


class PerfExperiment(Experiment):
    # Indicator specs are all prepened with "a_" to indicate they are after generation indicators
    INDICATOR_SPECS: ClassVar[list[tuple[str, bool]]] = [
        ("a_coverage", True),
        ("a_size", True),
        ("a_score", True),
        ("a_avg_dist_to_niches", False),
        ("a_avg_niche_count", True),
        ("a_filled_niches", True),
        ("a_avg_dist_to_ideal", False),
        ("a_niche_perf_avg", False),
    ]

    @staticmethod
    def parser(mode: ExperimentMode = ExperimentMode.RUN) -> argparse.ArgumentParser:
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
        parser.add_argument(
            "--kernel",
            action="store_true",
            help="Whether to use kernel-based pred group ornament",
        )
        parser.add_argument(
            "--abs",
            action="store_true",
            help="Whether to use absolute pred group ornament",
        )
        parser.add_argument(
            "--rel_fit",
            action="store_true",
            help="Whether to use relative fitness pred group ornament",
        )
        parser.add_argument(
            "--state",
            action="store_true",
            help="Whether to use state pred group ornament",
        )
        parser.add_argument(
            "--age", action="store_true", help="Whether to use age pred group ornament"
        )
        parser.add_argument(
            "--no_x",
            action="store_true",
            help="Whether to exclude the X feature from the ornament",
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
            "kernel",
            "abs",
            "rel_fit",
            "state",
            "age",
            "no_x",
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
        assert batch.problem is not None, "Batch problem is None"
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
        cols_top_keep = [ind for ind, _ in PerfExperiment.INDICATOR_SPECS] + [
            "a_fevals"
        ]
        batch_dfs = []
        pq_folder = Path(batches[0].parent_dir) / f"pq"
        pq_folder.mkdir(parents=True, exist_ok=True)
        out_files = []
        # Merge the record_stats files into one
        for batch in batches:
            file = PerfExperiment.file_paths(batch)["record_stats"]
            df = pd.read_csv(file).assign(seed=batch.seed)
            # only keep columns with indicators (and a_fevals)
            result_data = df[cols_top_keep]
            # group by a_fevals and take the mean of the indicators
            # This is because the data also contains x and y values, but the indicators are only computed after the generation, so we can just take the mean of the indicators for each generation
            result_data = result_data.groupby("a_fevals").mean().reset_index()
            result = CVResult(data=result_data, batch=batch, fevals_column="a_fevals")
            out_files.append(str(result.write(pq_folder)))
            batch_dfs.append(df)

        merged_record_stats = pd.concat(batch_dfs, ignore_index=True)

        # Save the merged record_stats to a new file
        out_path = Path(batches[0].parent_dir) / f"merged_record_stats_{exp_id}.csv"
        merged_record_stats.to_csv(out_path, index=False)
        out_files.append(str(out_path))
        return [
            out_files
        ]  # Return a list of lists of file paths, to be consistent with _run_batches

    @staticmethod
    def _produce_results(config: ResultsConfig, **kwargs) -> list[str]:
        # Build up result set from all input directories, i.e. across all runs of the experiment
        res_set = ResultSet()
        for in_dir in config.input_dirs:
            in_res_set = ResultSet.from_directory(in_dir)
            # rename the algorithm name in each result to the directory name, so that the plots can be distinguished by a human readable name
            algorithm_name = Path(in_dir).parent.name
            for res in in_res_set:
                res.algorithm = algorithm_name
            res_set.extend(in_res_set)
        plots = plot_rt(
            results=res_set,
            indicator_specs=PerfExperiment.INDICATOR_SPECS,
            output_dir=config.out_dir,
            file_prefix="mk_perf",
        )
        return plots

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
        problem: StaticProblem | None = None,
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
        elif pref_model == "ornament":
            from mating_kernel.pymoo.pref_model.ornament.dist_scorer import (
                DistScorer,
            )
            from mating_kernel.pymoo.pref_model.ornament.feature_parser import (
                FeatureOrnamentParser,
            )
            from mating_kernel.pymoo.pref_model.ornament_pref_model import (
                OrnamentPreferenceModel,
            )

            features = FeatureOrnamentParser.get_features(
                problem=problem,
                **pref_params,
            )

            pref_mod = OrnamentPreferenceModel(
                ornament_parser=FeatureOrnamentParser(features=features),
                scorer=DistScorer(),
            )
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
        problem: StaticProblem | None,
        algorithm_name: str,
        pref_model: str | None,
        pref_params: dict | None,
        num_candidates: int,
        candidate_pressure: int,
        n_gen: int,
        algorithm_params: dict | None = None,
        seed: int | None = None,
    ) -> tuple[Result, pd.DataFrame]:
        # NOTE: In the future, potentially identify problem via name and load it from a registry, rather than passing the problem object directly.
        if problem is None:
            raise ValueError("Problem must be provided for the experiment.")
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
            problem=problem,
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
