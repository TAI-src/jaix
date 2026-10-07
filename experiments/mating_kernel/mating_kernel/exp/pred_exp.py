import argparse
from mating_kernel.exp.exp import Experiment
import logging
from pathlib import Path
import pickle
import pandas as pd

from mating_kernel.exp.utils.batch import Batch
from mating_kernel.exp.utils.pre_pred_exp import (
    get_features,
    preprocess,
    get_data,
    target_types,
)
from mating_kernel.exp.utils.pred import run_analysis

logger = logging.getLogger(__name__)


# TODO: Add information computed on ideal point later


class PredExperiment(Experiment):
    @staticmethod
    def parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument(
            "--perf_stats_dir",
            type=str,
            default=None,
            help="Directory to find the performance statistics for the prediction experiment.",
        )
        parser.add_argument(
            "--kernel",
            action="store_true",
            help="Wether to use the kernel information for the prediction experiment.",
        )
        parser.add_argument(
            "--abs",
            action="store_true",
            help="Whether to use the absolute information for the prediction experiment.",
        )

        parser.add_argument(
            "--rel_fit",
            action="store_true",
            help="Whether to use the fitness information for the prediction experiment.",
        )
        parser.add_argument(
            "--state",
            action="store_true",
            help="Whether to use the state information of the populations",
        )
        parser.add_argument(
            "--age",
            action="store_true",
            help="Whether to use the age information for the prediction experiment.",
        )
        parser.add_argument(
            "--keep_mutated",
            action="store_true",
            help="Whether to keep mutated children in the prediction experiment.",
        )
        parser.add_argument(
            "--no_x",
            action="store_true",
            help="Whether to exclude the x information for the prediction experiment.",
        )

        parser.add_argument(
            "--target",
            type=str,
            nargs="+",
            choices=["survived", "o_dist_to_ideal", "o_F_0", "o_F_1"],
            default="survived",
            help="The target variable for the prediction experiment.",
        )
        parser.add_argument(
            "--feature_analysis",
            action="store_true",
            help="Whether to perform feature analysis for the prediction experiment.",
        )

        return parser

    # TODO: This should somehow move up to general Experiment
    @staticmethod
    def _run_batches(batches: list[Batch], **kwargs) -> list[list[str]]:
        batch_files = []
        for b in batches:
            logger.debug(f"Running batch {b.name}")
            files = PredExperiment._run_batch(b, **kwargs)
            batch_files.append(files)
        return batch_files

    @staticmethod
    def file_paths(batch: Batch) -> dict[str, Path]:
        file_paths_dict = {
            "dataset": batch.out_dir / f"dataset_{batch.name}.csv",
            "result": batch.out_dir / f"result_{batch.name}.pkl",
        }
        return file_paths_dict

    @staticmethod
    def setting_name(batch: Batch) -> str:
        """
        Return a string representation of the batch settings by giving back only true settings plus target
        """
        settings = []
        if batch.kernel:  # type: ignore[attr-defined]
            settings.append("kernel")
        if batch.abs:  # type: ignore[attr-defined]
            settings.append("abs")
        if batch.rel_fit:  # type: ignore[attr-defined]
            settings.append("rel_fit")
        if batch.state:  # type: ignore[attr-defined]
            settings.append("state")
        if batch.age:  # type: ignore[attr-defined]
            settings.append("age")
        if batch.keep_mutated:  # type: ignore[attr-defined]
            settings.append("keep_mutated")
        if batch.no_x:  # type: ignore[attr-defined]
            settings.append("no_x")
        settings.append(batch.target)  # type: ignore[attr-defined]
        return "_".join(settings)

    @staticmethod
    def _post_process_batches(batches: list[Batch], **kwargs) -> list[list[str]]:
        # Create dataframe with all cv results and save to csv
        result_files = [PredExperiment.file_paths(batch)["result"] for batch in batches]
        names = [PredExperiment.setting_name(batch) for batch in batches]
        # Load all results
        results = {}
        for name, result_file, batch in zip(names, result_files, batches):
            with open(result_file, "rb") as f:
                res = pickle.load(f)
            results[(name, batch.pid)] = {
                "cv_score_mean": res.cv_score_mean,
                "cv_score_std": res.cv_score_std,
            }
        df = pd.DataFrame.from_dict(results, orient="index")
        df.index = pd.MultiIndex.from_tuples(df.index, names=["setting", "pid"])
        df = df.reset_index()
        out_dir = batches[0].out_dir.parent
        res_names = "_".join(sorted(set(names)))
        out_file = out_dir / f"summary_{res_names}.csv"
        df.to_csv(out_file)

        # combine all summary files in the same directory into one summary file
        # FIXME: A bit dirty, but needs this until we can combine true/false batches
        summary_file = out_dir / "cv_results_agg.csv"
        if summary_file.exists():
            df.to_csv(summary_file, mode="a", header=False, index=False)
        else:
            df.to_csv(summary_file, index=False)

        return [[str(out_file), str(summary_file)]]

    @staticmethod
    def _run_batch(batch: Batch, **kwargs) -> list[str]:
        out_files = PredExperiment.file_paths(batch)
        # Find the file for the batch
        df = get_data(batch.perf_stats_dir, batch.pid)  # type: ignore[attr-defined]
        # Determine the features to use
        features = get_features(
            batch.kernel,  # type: ignore[attr-defined]
            batch.abs,  # type: ignore[attr-defined]
            batch.rel_fit,  # type: ignore[attr-defined]
            batch.state,  # type: ignore[attr-defined]
            batch.age,  # type: ignore[attr-defined]
            batch.no_x,  # type: ignore[attr-defined]
        )
        # Preprocess
        df = preprocess(
            features=features + [batch.target],  # type: ignore[attr-defined]
            df=df,
            remove_mutated=not batch.keep_mutated,  # type: ignore[attr-defined]
        )
        df.to_csv(out_files["dataset"], index=False)
        # run the analysis
        res = run_analysis(
            df,
            input_cols=list(df.columns.difference([batch.target])),  # type: ignore[attr-defined]
            target_col=batch.target,  # type: ignore[attr-defined]
            target_type=target_types[batch.target],  # type: ignore[attr-defined]
            group_cols=["seed"],
            skip_feature_analysis=not batch.feature_analysis,  # type: ignore[attr-defined]
            random_state=batch.seed,
            **kwargs,
        )
        # Save the result
        with open(out_files["result"], "wb") as f:
            pickle.dump(res, f)

        return [str(f) for f in out_files.values()]
