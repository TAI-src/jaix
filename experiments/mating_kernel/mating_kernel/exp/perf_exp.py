from pymoo.algorithms.moo.nsga2 import NSGA2

from pymoo.optimize import minimize
from mating_kernel.pymoo.parser.population_parser import PopulationParser
from mating_kernel.pymoo.parser.reproduction_parser import ReproductionParser
from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper
from mating_kernel.pymoo.recordable_object import make_recordable
from mating_kernel.pymoo.recording_callback import RecordingCallback
from mating_kernel.experiments.experiment import Experiment, ExperimentConfig
from mating_kernel.exp.utils.parse_args import SettingArg
from mating_kernel.exp.utils.batch import Batch
from jaix.env.utils.problem.static_problem import StaticProblem
from typing import Any
import argparse


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
            algorithm_name=batch.algorithm_name,
            selector=batch.selector,
            n_gen=batch.n_gen,
            algorithm_params={},
            seed=batch.seed,
        )
        # Save the result and record_stats to files in batch.out_dir
        # TODO:what to save here?
        return []

    @staticmethod
    def run_instrumented_pymoo(
        problem: StaticProblem,
        algorithm_name: str,
        selector: str,
        n_gen: int,
        algorithm_params: dict = {},
        seed: int | None = None,
    ) -> tuple[Any, list[dict]]:
        parser_classes = [ReproductionParser, PopulationParser]
        record_args = list(set(p.record_args for p in parser_classes))
        record_attributes = list(set(p.record_attributes for p in parser_classes))
        if algorithm_name == "NSGA2":
            algorithm_class = NSGA2
        else:
            raise ValueError(f"Unsupported algorithm: {algorithm_name}")
        recorded_algorithm_class = make_recordable(algorithm_class)
        if selector is not None:
            raise NotImplementedError("Selector conversion is not implemented yet.")
        algorithm = recorded_algorithm_class(
            **algorithm_params,
            record_args=record_args,
            record_attributes=record_attributes,
        )
        pymoo_problem = PymooProblemWrapper(problem)
        callback = RecordingCallback(
            recording_parsers=[p(ideal=problem.ideal_point) for p in parser_classes]
        )
        result = minimize(
            pymoo_problem,
            algorithm,
            seed=seed,
            termination=("n_gen", n_gen),
            callback=callback,
        )
        return result, callback.data["record_stats"]
