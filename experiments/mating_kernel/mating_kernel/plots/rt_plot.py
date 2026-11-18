import logging
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from cocoviz import (
    Indicator,
    ProblemDescription,
    ResultSet,
    rtpplot,
)


from mating_kernel.problems.problem_info import ProblemInfo
from mating_kernel.exp.utils.batch import Batch
from mating_kernel.plots.cv_result import Result

logger = logging.getLogger(__name__)


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
INDICATOR_SPECS = [(f"a_{ind}", lib) for ind, lib in indicator_specs]


def result_from_df(
    source_df: pd.DataFrame,
    batch: Batch,
    indicators: list[str] | None,
) -> Result:
    if indicators is None:
        indicators = [ind for ind, _ in INDICATOR_SPECS]
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
    # Check that we only have one seed
    assert source_df["seed"].nunique() == 1

    # group dataset by seed and a_fevals, and compute the mean of each indicator for each group
    indicator_means = source_df.groupby(["a_fevals"])[indicators].mean().reset_index()
    result = Result(
        data=indicator_means,
        batch=batch,
        fevals_column="a_fevals",
    )
    return result


def get_result_set(
    source_df: pd.DataFrame | None,
    indicators: list[str],
    pq_path: Path | None = None,
) -> ResultSet:
    # Source df is not provided, try reading existing results
    if source_df is None:
        if pq_path is None:
            raise ValueError("source_df is None and pq_path is None. Cannot load data.")
        # Check if we can load from the pq_path
        logger.info(f"source_df is None, attempting to load from pq_path {pq_path}")
        res_set = ResultSet.from_directory(pq_path)
        return res_set
    # Source df is provided, get results from source_df
    assert "seed" in source_df.columns and "exp_id" in source_df.columns
    res_set = ResultSet()
    for df in source_df.groupby(["exp_id", "seed"]):
        res = result_from_df(source_df, batch, indicators)
        res_set.append(res)

    return res_set


def write_result_set(
    res_set: ResultSet,
    pq_path: Path | str,
) -> list[str]:
    file_names = []
    for result in res_set:
        file_name = str(pq_path / result.pq_file_name)
        result.to_parquet(file_name)
    return file_names


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
        indicator_specs = INDICATOR_SPECS
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
