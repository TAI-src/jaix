import argparse
import copy
import logging
import pickle

import pandas as pd
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.algorithm import Algorithm
from pymoo.core.result import Result
from pymoo.optimize import minimize

from jaix.env.utils.problem.static_problem import StaticProblem
from mating_kernel.exp.exp import Experiment
from mating_kernel.exp.utils.batch import Batch
from mating_kernel.pymoo.offspring_success_recording_callback import (
    OffspringSuccessRecordingCallback,
)
from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper
from mating_kernel.pymoo.recordable_object import make_recordable

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
            "--selector",
            nargs="+",
            type=str,
            required=False,
            default=None,
        )
        parser.add_argument(
            "--n_gen",
            type=int,
            default=1000,
            help="Number of generations to run the algorithm for.",
        )
        return parser

    @staticmethod
    def _run_batch(batch: Batch, **kwargs) -> list[str]:
        result, record_stats = PerfExperiment.run_instrumented_pymoo(
            problem=batch.problem,
            algorithm_name=batch.alg_name,  # type: ignore[attr-defined]
            selector=batch.selector,  # type: ignore[attr-defined]
            n_gen=batch.n_gen,  # type: ignore[attr-defined]
            algorithm_params={},  # type: ignore[attr-defined]
            seed=batch.seed,
        )
        # Save the result and record_stats to files in batch.out_dir
        result_file = batch.out_dir / f"result_{batch.name}.pkl"
        with open(result_file, "wb") as f:
            result_cpy = copy.deepcopy(result)
            result_cpy.problem = str(batch.problem)
            result_cpy.algorithm = batch.alg_name  # type: ignore[attr-defined]
            pickle.dump(result_cpy, f)
        record_stats_file = batch.out_dir / f"record_stats_{batch.name}.csv"
        record_stats.to_csv(record_stats_file, index=False)
        batch_file = batch.out_dir / f"batch_{batch.name}.pkl"
        with open(batch_file, "wb") as f:
            batch_cpy = batch.model_copy()
            batch_cpy.problem = str(batch.problem)
            pickle.dump(batch_cpy, f)

        return [str(result_file), str(record_stats_file), str(batch_file)]

    @staticmethod
    def get_alg_class(algorithm_name: str) -> Algorithm:
        if algorithm_name == "NSGA2":
            return NSGA2
        else:
            raise ValueError(f"Unsupported algorithm: {algorithm_name}")

    @staticmethod
    def get_recorded_alg(
        algorithm_name: str,
        selector: str | None,
        algorithm_params: dict | None,
        record_args: list[str] | None = None,
        record_attributes: list[str] | None = None,
    ) -> Algorithm:
        algorithm_class = PerfExperiment.get_alg_class(algorithm_name)
        record_alg_class = make_recordable(algorithm_class)
        if algorithm_params is None:
            algorithm_params = {}
        if selector is not None:
            raise NotImplementedError("Selector conversion is not implemented yet.")
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
        selector: str | None,
        n_gen: int,
        algorithm_params: dict | None = None,
        seed: int | None = None,
    ) -> tuple[Result, pd.DataFrame]:
        pymoo_problem = PymooProblemWrapper(problem)
        callback = OffspringSuccessRecordingCallback(pymoo_problem)
        algorithm = PerfExperiment.get_recorded_alg(
            algorithm_name=algorithm_name,
            selector=selector,
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
