from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from cocoviz import (
    Indicator,
    ProblemDescription,
    Result as CVResult,
    ResultSet,
    rtpplot,
)
from cocoviz._typing import FrameLike

from mating_kernel.problems.problem_info import ProblemInfo
import logging

logger = logging.getLogger(__name__)


class Result(CVResult):
    def __init__(
        self,
        algorithm: str,
        problem: ProblemDescription,
        data: FrameLike,
        fevals_column: str = "fevals",
        seed: int | None = None,
    ):
        self.seed = seed
        self.pq_file_name = f"{problem.name}_{algorithm}_{seed}.parquet"
        super().__init__(algorithm, problem, data, fevals_column)


def get_description(pinfo: ProblemInfo) -> ProblemDescription:
    cocoviz_description = ProblemDescription(
        name=pinfo.uuid,
        instance=0,
        number_of_variables=pinfo.num_variables,
        number_of_objectives=pinfo.num_objectives,
    )
    return cocoviz_description


def results_from_df(
    source_df: pd.DataFrame,
    exp_name: str,
    pdesc: ProblemDescription,
    indicators: list[str],
) -> ResultSet:

    # Check that all indicators are present in the source_df and are numeric
    assert all(
        indicator in source_df.columns
        and pd.api.types.is_numeric_dtype(source_df[indicator])
        for indicator in indicators
    ), "All indicators must be present in the source_df and must be numeric."

    # Check that seed and a_fevals columns are present in the source_df
    assert (
        "seed" in source_df.columns and "a_fevals" in source_df.columns
    ), "source_df must contain 'seed' and 'a_fevals' columns."

    # group dataset by seed and a_fevals, and compute the mean of each indicator for each group
    indicator_means = (
        source_df.groupby(["seed", "a_fevals"])[indicators].mean().reset_index()
    )
    result_set = ResultSet()
    for seed, run_data in indicator_means.groupby("seed"):
        result = Result(
            algorithm=exp_name,
            problem=pdesc,
            data=run_data,
            fevals_column="a_fevals",
            seed=seed,
        )
        result_set.append(result)
    return result_set


def get_result_set(
    source_df: pd.DataFrame | None,
    exp_name: str,
    pdesc: ProblemDescription,
    indicators: list[str],
    pq_path: Path | None = None,
) -> ResultSet:
    if source_df is None:
        if pq_path is None:
            raise ValueError("source_df is None and pq_path is None. Cannot load data.")
        # Check if we can load from the pq_path
        logger.info(f"source_df is None, attempting to load from pq_path {pq_path}")
        res_set = ResultSet.from_directory(pq_path)
        return res_set
    # Source df is provided, update results in pq_path if they exist
    res_set = results_from_df(source_df, exp_name, pdesc, indicators)

    if pq_path is not None:
        # Update the results in the pq_path if they
        logger.info(
            f"Writing results to pq_path {pq_path}. This will overwrite existing files with the same name."
        )
        for result in res_set:
            result.to_parquet(pq_path / result.pq_file_name)
    return res_set


def plot_rt(
    results: ResultSet,
    indicator_specs: list[tuple[str, bool]],
    output_dir: str | Path,
    file_prefix: str = "",
) -> list[str]:
    file_names = []
    for problem, by_prob in results.by_problem_name():
        for ind_col, lib in indicator_specs:
            ind = Indicator(ind_col, larger_is_better=lib)
            ax = rtpplot(results, ind)
            # save the plot to a file
            file_name = (
                f"{file_prefix}_rt_{problem}_{ind_col}.pdf"
                if file_prefix
                else f"rt_{problem}_{ind_col}.pdf"
            )
            fig = ax.figure
            file_path = Path(output_dir) / Path(file_name)
            fig.savefig(file_path, dpi=300)
            plt.close(fig)
            file_names.append(str(file_path))
    return file_names


def plot_rt_cocoviz(
    source: pd.DataFrame | Path | str,
    exp_name: str,
    pinfo: ProblemInfo,
    output_dir: str | Path,
    indicator_specs: list[tuple[str, bool]] | None = None,
    file_prefix: str = "",
) -> list[str]:
    if indicator_specs is None:
        indicator_specs = [
            ("coverage", True),
            ("size", True),
            ("score", True),
            ("avg_dist_to_niches", False),
            ("avg_niche_count", True),
            ("filled_niches", True),
            ("avg_dist_to_ideal", False),
            ("niche_perf_avg", False),
        ]
        # prepend "a_" to each indicator name (for after the generation)
        indicator_specs = [(f"a_{ind}", lib) for ind, lib in indicator_specs]
    pdesc = get_description(pinfo)
    if isinstance(source, pd.DataFrame):
        # load result set from source_df
        result_set = get_result_set(
            source_df=source,
            exp_name=exp_name,
            pdesc=pdesc,
            indicators=[ind for ind, _ in indicator_specs],
            pq_path=Path(output_dir),
        )
    else:
        # Can read from source path
        assert Path(source).exists()
        result_set = ResultSet.from_directory(source)

    file_names = plot_rt(result_set, indicator_specs, output_dir, file_prefix)
    return file_names
