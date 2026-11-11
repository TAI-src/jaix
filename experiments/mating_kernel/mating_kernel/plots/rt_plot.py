from cocoviz import ProblemDescription, Result, ResultSet, Indicator, rtpplot
from mating_kernel.problems.problem_info import ProblemInfo
import pandas as pd

from pathlib import Path
import matplotlib.pyplot as plt


def get_description(pinfo: ProblemInfo) -> ProblemDescription:
    cocoviz_description = ProblemDescription(
        name=pinfo.uuid,
        instance=0,
        number_of_variables=pinfo.num_variables,
        number_of_objectives=pinfo.num_objectives,
    )
    return cocoviz_description


def get_result_set(
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
            fevals_col="a_fevals",
        )
        result_set.append(result)
    return result_set


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
            fig.savefig(file_name, dpi=300)
            plt.close(fig)
            file_names.append(file_name)
    return file_names


def plot_rt_cocoviz(
    source_df: pd.DataFrame,
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
    result_set = get_result_set(
        source_df, exp_name, pdesc, [ind for ind, _ in indicator_specs]
    )
    file_names = plot_rt(result_set, indicator_specs, output_dir, file_prefix)
    return file_names
