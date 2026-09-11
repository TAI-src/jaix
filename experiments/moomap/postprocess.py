import os
from pathlib import Path
from collections import defaultdict
from nsga3_experiment import NSGA3ExperimentConfig
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from cobi_config_generator import old_names

import json


def build_results_dict(results_dir: str | None = None) -> dict:
    """
    Build a dictionary of results from the logs directory.
    """
    if results_dir is None:
        results_dir = f"{os.path.dirname(os.path.abspath(__file__))}/results"
    res_dict = defaultdict(dict)
    # Get all csv files in the logs directory
    csv_files = list(Path(results_dir).rglob("*.csv"))
    for file in csv_files:
        problem_str = file.name.replace("results_", "").replace(".csv", "")
        tmp_dict = {"file": file}

        # get example config to know number of niches (dirty)
        # remove results_ prefix and .csv suffix from problem name to get config name
        config_file = file.parent / f"config_{problem_str}.json"
        # load the config to get the number of niches

        with open(config_file, "r") as f:
            config = json.load(f)
            tmp_dict["config"] = config["NSGA3ExperimentConfig"]

        res_dict[problem_str][file.parent.name] = tmp_dict

    res_dict = dict(res_dict)  # convert to normal dict
    expected_problem_list = NSGA3ExperimentConfig.generate_problem_list()
    problem_strs = [str(problem) for problem in expected_problem_list]
    for problem_str, problem in zip(problem_strs, expected_problem_list):
        if problem_str.startswith("cobi_"):
            # translate old names to new names
            problem_str = old_names.get(problem_str, problem_str)
        if problem_str not in res_dict:
            print(f"Missing results for problem: {problem_str}")
        else:
            # add info about the problem to the dict
            res_dict[problem_str]["ideal_point"] = problem.ideal_point
            res_dict[problem_str]["nadir_point"] = problem.nadir_point
            res_dict[problem_str]["num_objectives"] = problem.num_objectives
            res_dict[problem_str]["num_variables"] = problem.dimension
            res_dict[problem_str]["lower_bounds"] = problem.lower_bounds
            res_dict[problem_str]["upper_bounds"] = problem.upper_bounds

    return res_dict


def get_success_cols(df: pd.DataFrame) -> list:
    """
    Get the columns related to success rate of offspring.
    Adds additional normalised versions of offspring rank as well based on generation
    """
    # normalise offspring rank by generation
    # rank 0 stays 0, worst rank is 1, and all other ranks are scaled between 0 and 1
    df["n_offspring_rank"] = df.groupby(["generation", "seed"])[
        "offspring_rank"
    ].transform(lambda x: x / x.max())
    df.loc[df["offspring_rank"] == 0, "n_offspring_rank"] = 0
    # critical rank is the worst (i.e. max) rank achieved in the generation, which is stored in archive_stats_after_max_rank
    df["critical_rank"] = df["archive_stats_after_max_rank"]

    # normalise offspring rank by critical rank per generation
    df["ncrit_offspring_rank"] = df["offspring_rank"] / df["critical_rank"]
    # if critical rank is 0, set ncrit_offspring_rank to 0
    df.loc[df["critical_rank"] == 0, "ncrit_offspring_rank"] = 0

    success_cols = [
        "offspring_added",
        "offspring_rank",
        "n_offspring_rank",
        "ncrit_offspring_rank",
        "offspring_dist_to_ideal",
    ]
    return success_cols


def compile_niche_success(
    data: pd.DataFrame, num_parents: int = 2, per_niche: bool = False
) -> pd.DataFrame:
    """
    Compile the success rate of each niche pair in the data.
    """
    df = data.copy()
    # get different pairings of parent niches and their success rate
    parent_niche_cols = ["parent_{i}_niche".format(i=i) for i in range(num_parents)]
    # sort the parent niche columns to ensure that the order of parents does not matter
    df[parent_niche_cols] = np.sort(df[parent_niche_cols].values, axis=1)
    if not per_niche:
        group_cols = parent_niche_cols
    else:
        group_cols = parent_niche_cols + ["offspring_niche"]
    success_cols = get_success_cols(df)
    niche_success = (
        df.groupby(group_cols)[success_cols]
        .agg(
            offspring_added_rate=("offspring_added", "mean"),
            offspring_rank_mean=("offspring_rank", "mean"),
            offspring_rank_std=("offspring_rank", "std"),
            num_offspring=("offspring_added", "count"),
            n_offspring_rank_mean=("n_offspring_rank", "mean"),
            n_offspring_rank_std=("n_offspring_rank", "std"),
            ncrit_offspring_rank_mean=("ncrit_offspring_rank", "mean"),
            ncrit_offspring_rank_std=("ncrit_offspring_rank", "std"),
            offspring_dist_to_ideal_mean=(
                "offspring_dist_to_ideal",
                "mean",
            ),  # This auto-removes nans
        )
        .reset_index()
    )

    return niche_success


def plot_niche_success(
    niche_success: pd.DataFrame,
    max_niche: int,
    output_dir: str,
    success_col: str = "offspring_added_rate",
    file_prefix: str = "",
) -> None:
    """
    Plot the success rate of each niche pair in the data.
    """

    # create a heatmap of the offspring_added_rate for each niche pair
    # niches should be from 0 to max_niche, so we can create a pivot table with all combinations of niches
    idx = pd.MultiIndex.from_product(
        [range(max_niche), range(max_niche)], names=["parent_0_niche", "parent_1_niche"]
    )
    niche_success = (
        niche_success.set_index(["parent_0_niche", "parent_1_niche"])
        .reindex(idx)
        .reset_index()
    )
    pivot_table = niche_success.pivot(
        index="parent_0_niche", columns="parent_1_niche", values=success_col
    )
    cmap = sns.color_palette("viridis", as_cmap=True)
    cmap.set_bad("lightgray")
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        pivot_table,
        annot=False,
        fmt=".2f",
        cmap=cmap,
        cbar_kws={"label": success_col},
    )
    plt.gca().invert_yaxis()
    plt.title(f"{file_prefix}: {success_col} by Parent Niche Pair")
    plt.xlabel("Parent 1 Niche")
    plt.ylabel("Parent 0 Niche")
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, f"{file_prefix}niche_success_{success_col}.png")
    )
    plt.close()


def angle_between(
    vec1: np.ndarray, vec2: np.ndarray, axis: int = -1, min_norm: float = 1e-12
) -> np.ndarray | np.floating:
    """
    Compute the angle between two vectors in radians.
    """
    norm1 = np.linalg.norm(vec1, axis=axis)
    norm2 = np.linalg.norm(vec2, axis=axis)

    valid = (norm1 > min_norm) & (norm2 > min_norm)

    dot_product = np.sum(vec1 * vec2, axis=axis)
    cos_angle = np.divide(
        dot_product,
        norm1 * norm2,
        out=np.full_like(dot_product, np.nan, dtype=float),
        where=valid,
    )

    return np.arccos(np.clip(cos_angle, -1.0, 1.0))


def compile_distance_success(
    data: pd.DataFrame, num_parents: int = 2, ideal_point: np.ndarray | None = None
) -> pd.DataFrame:
    """
    Get data related to success rate of offspring based on the distance between parents.
    """
    assert num_parents == 2, "This function only supports 2 parents."
    df = data.copy()
    # compute parent distances in x and y
    for dim in ["x", "y"]:
        parent_cols = [
            "parent_{i}_{dim}".format(i=i, dim=dim) for i in range(num_parents)
        ]
        for col in parent_cols:
            df[col] = df[col].apply(lambda x: np.fromstring(x.strip("[]"), sep=" "))
        df["parent_{dim}_distance".format(dim=dim)] = np.linalg.norm(
            np.stack(df[parent_cols[0]].values) - np.stack(df[parent_cols[1]].values),
            axis=1,
        )

    success_cols = get_success_cols(df)
    additional_cols = [
        "generation",
        "parent_x_distance",
        "parent_y_distance",
        "parent_0_dist_to_ideal",
        "parent_1_dist_to_ideal",
    ]
    # add angle between vectors of parents on y dim with ideal point as origin
    if ideal_point is not None:
        vec1 = np.stack(df["parent_0_y"].values) - ideal_point
        vec2 = np.stack(df["parent_1_y"].values) - ideal_point
        df["parent_angle"] = angle_between(vec1, vec2)
        additional_cols.append("parent_angle")

    # return the relevant columns for distance success analysis
    distance_success = df[additional_cols + success_cols]
    return distance_success


def plot_distance_success_2d(
    distance_success: pd.DataFrame,
    output_dir: str,
    x_axis="parent_x_distance",
    y_axis="parent_y_distance",
    success_col: str = "offspring_dist_to_ideal",
    file_prefix: str = "",
) -> None:
    """
    Plot the success rate of offspring based on the distance between parents.
    """
    # filter rows with nan values in success_col
    # This happens because some scores are not computed for offspring if they are not added

    dist_succ = distance_success.copy()
    dist_succ = dist_succ.dropna(subset=[success_col])

    norm = mpl.colors.Normalize(
        vmin=min(dist_succ[success_col]), vmax=max(dist_succ[success_col])
    )
    cmap = plt.get_cmap("viridis")

    plt.figure(figsize=(10, 8))
    ax = sns.scatterplot(
        data=dist_succ,
        x=x_axis,
        y=y_axis,
        hue=success_col,
        palette=cmap,
        hue_norm=norm,
        alpha=0.7,
        legend=False,
    )
    # Create continuous colorbar
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    plt.colorbar(sm, ax=ax, label=success_col)

    only_added = len(distance_success) > len(dist_succ)
    if only_added:
        plt.title(f"{file_prefix}: {success_col} by {x_axis} vs {y_axis} (only added)")
    else:
        plt.title(f"{file_prefix}: {success_col} by {x_axis} vs {y_axis}")
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    # plt.colorbar(label=success_col)
    plt.tight_layout()
    plt.savefig(
        os.path.join(
            output_dir, f"{file_prefix}dis_{x_axis}vs{y_axis}_{success_col}.png"
        )
    )
    plt.close()


def plot_distance_success_1d(
    distance_success: pd.DataFrame,
    output_dir: str,
    dist_col: str = "parent_x_distance",
    success_col: str = "offspring_dist_to_ideal",
    file_prefix: str = "",
) -> None:
    """
    Plot the success rate of offspring based on the distance between parents in 1D.
    Use the generation as the colour of the points to see if there is a trend over time.
    """
    norm = mpl.colors.Normalize(
        vmin=min(distance_success["generation"]),
        vmax=max(distance_success["generation"]),
    )
    cmap = plt.get_cmap("viridis")
    plt.figure(figsize=(10, 8))
    ax = sns.scatterplot(
        data=distance_success,
        x=dist_col,
        y=success_col,
        hue="generation",
        palette=cmap,
        hue_norm=norm,
        alpha=0.7,
        legend=False,
    )

    # Create continuous colorbar
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    plt.colorbar(sm, ax=ax, label="Generation")

    plt.title(f"{file_prefix}: {success_col} by Parent Distance ({dist_col})")
    plt.xlabel(f"Parent Distance ({dist_col})")
    plt.ylabel(success_col)
    plt.tight_layout()
    plt.savefig(
        os.path.join(
            output_dir,
            f"{file_prefix}_gen_distance_success_{dist_col}_{success_col}.png",
        )
    )
    plt.close()


def plot_sankey_flows(
    niche_success: pd.DataFrame,
    output_dir: str,
    file_prefix: str = "",
    num_top_parent_flows: int = 20,
    num_top_child_flows: int = 5,
    max_width: float = 0.025,
    min_width: float = 0.003,
    node_height: float = 0.035,
    node_width: float = 0.025,
    alpha: float = 0.5,
    left_x=0.05,
    right_x=0.95,
    bezier_control_factor: float = 0.4,
    flow_colour: str = "steelblue",
    source_colour: str = "black",
    target_colour: str = "darkorange",
    legend_x=0.35,
    legend_y=0.04,
    legend_spacing=0.035,
) -> None:
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch

    # preprocess df
    df = niche_success.copy()
    sel_cols = ["parent_0_niche", "parent_1_niche", "offspring_niche", "num_offspring"]
    df = df[sel_cols]
    # remove columns where offspring_niche == -1 (i.e. no offspring added)
    df = df[df["offspring_niche"] != -1]
    # Add parent pair as a new column
    df["parent_pair"] = df.apply(
        lambda row: f"{row['parent_0_niche']}-{row['parent_1_niche']}", axis=1
    )

    # Aggregate flows
    df_grouped = (
        df.groupby(["parent_pair", "offspring_niche"])["num_offspring"]
        .sum()
        .reset_index()
    )
    # Total offspring per parent pair
    parent_totals = (
        df_grouped.groupby("parent_pair")["num_offspring"]
        .sum()
        .sort_values(ascending=False)
    )

    # Select the most important parent pairs
    top_parent_pairs = parent_totals.head(num_top_parent_flows).index

    df_selected = df_grouped[df_grouped["parent_pair"].isin(top_parent_pairs)].copy()

    # For each parent pair, keep its top K offspring flows
    df_selected = (
        df_selected.sort_values(
            ["parent_pair", "num_offspring"],
            ascending=[True, False],
        )
        .groupby("parent_pair")
        .head(num_top_child_flows)
    )
    df_grouped = df_selected.copy()

    # ---------------------------------------------------------
    # Nodes
    # ---------------------------------------------------------
    parent_pairs = df_grouped["parent_pair"].unique().tolist()
    offspring_niches = sorted(df_grouped["offspring_niche"].unique())

    # Assign vertical positions
    # Top = 1, bottom = 0
    parent_y = {
        pair: 1 - (i + 1) / (len(parent_pairs) + 1)
        for i, pair in enumerate(parent_pairs)
    }

    offspring_y = {
        niche: 1 - (i + 1) / (len(offspring_niches) + 1)
        for i, niche in enumerate(offspring_niches)
    }

    # ---------------------------------------------------------
    # Scale flow widths
    # ---------------------------------------------------------
    max_flow = df_grouped["num_offspring"].max()

    # Maximum vertical thickness of a flow

    # ---------------------------------------------------------
    # Figure
    # ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(14, 10))

    # ---------------------------------------------------------
    # Draw flows
    # ---------------------------------------------------------
    for _, row in df_grouped.iterrows():
        source = row["parent_pair"]
        target = row["offspring_niche"]
        value = row["num_offspring"]

        y0 = parent_y[source]
        y1 = offspring_y[target]

        width = max(min_width, max_width * value / max_flow)

        # Cubic Bezier control points
        dx = right_x - left_x

        control_x1 = left_x + bezier_control_factor * dx
        control_x2 = right_x - bezier_control_factor * dx

        # Upper boundary
        verts_top = [
            (left_x, y0 + width / 2),
            (control_x1, y0 + width / 2),
            (control_x2, y1 + width / 2),
            (right_x, y1 + width / 2),
        ]

        # Lower boundary
        verts_bottom = [
            (right_x, y1 - width / 2),
            (control_x2, y1 - width / 2),
            (control_x1, y0 - width / 2),
            (left_x, y0 - width / 2),
        ]

        vertices = verts_top + verts_bottom + [(left_x, y0 + width / 2)]

        codes = (
            [Path.MOVETO]
            + [Path.CURVE4] * 3
            + [Path.LINETO]
            + [Path.CURVE4] * 3
            + [Path.CLOSEPOLY]
        )

        path = Path(vertices, codes)

        patch = PathPatch(
            path,
            facecolor=flow_colour,
            edgecolor="none",
            alpha=alpha,
        )

        ax.add_patch(patch)

    # ---------------------------------------------------------
    # Flow scale legend
    # ---------------------------------------------------------
    min_flow = df_grouped["num_offspring"].min()
    max_flow = df_grouped["num_offspring"].max()
    mid_flow = (min_flow + max_flow) / 2

    legend_values = [min_flow, mid_flow, max_flow]

    ax.text(
        legend_x,
        legend_y + 0.10,
        "Number of offspring",
        fontsize=10,
        fontweight="bold",
    )

    for i, value in enumerate(legend_values):
        width = max_width * value / max_flow
        width = max(min_width, width)

        y = legend_y + (2 - i) * legend_spacing

        ax.plot(
            [legend_x, legend_x + 0.15],
            [y, y],
            color=flow_colour,
            linewidth=width * 400,
            alpha=alpha,
            solid_capstyle="butt",
        )

        ax.text(
            legend_x + 0.17,
            y,
            f"{value:,.0f}",
            va="center",
            fontsize=9,
        )

    # ---------------------------------------------------------
    # Draw nodes
    # ---------------------------------------------------------

    for pair, y in parent_y.items():
        ax.add_patch(
            plt.Rectangle(
                (left_x - node_width / 2, y - node_height / 2),
                node_width,
                node_height,
                color=source_colour,
            )
        )

        ax.text(
            left_x - node_width,
            y,
            pair,
            ha="right",
            va="center",
            fontsize=9,
        )

    for niche, y in offspring_y.items():
        ax.add_patch(
            plt.Rectangle(
                (right_x - node_width / 2, y - node_height / 2),
                node_width,
                node_height,
                color=target_colour,
            )
        )

        ax.text(
            right_x + 0.025,
            y,
            str(niche),
            ha="left",
            va="center",
            fontsize=9,
        )

    # ---------------------------------------------------------
    # Formatting
    # ---------------------------------------------------------
    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(0, 1)

    ax.text(
        left_x,
        1.03,
        "Parent pair",
        ha="center",
        va="bottom",
        fontsize=13,
        fontweight="bold",
    )

    ax.text(
        right_x,
        1.03,
        "Offspring niche",
        ha="center",
        va="bottom",
        fontsize=13,
        fontweight="bold",
    )

    ax.set_title(
        f"{file_prefix}: Top {num_top_parent_flows} Parent Pairs and Top {num_top_child_flows} Offspring Niches",
        fontsize=15,
    )

    ax.axis("off")

    plt.tight_layout()

    fig.savefig(
        os.path.join(
            output_dir,
            f"{file_prefix}_sankey_flows.pdf",
        ),
        bbox_inches="tight",
    )

    plt.close(fig)


def postprocess_results(problem_dict: dict, results_dir: str = ".") -> None:

    for problem, problem_dict in problem_dict.items():
        print(f"Processing results for problem: {problem}")

        runs = [key for key in problem_dict.keys() if key.startswith("r_")]
        res_files = [problem_dict[run]["file"] for run in runs]
        # FIXME: this is a workaround because we know that the number of niches is the same for all runs of a problem, but we should probably check that
        num_niches = problem_dict[runs[0]]["config"]["num_offspring"]

        # combine all the result files into a single dataframe
        problem_df = pd.concat([pd.read_csv(f) for f in res_files], ignore_index=True)
        per_niche_success = compile_niche_success(
            problem_df, num_parents=2, per_niche=True
        )
        plot_sankey_flows(
            per_niche_success,
            output_dir=results_dir,
            file_prefix=problem,
        )
        per_niche_success.to_csv(
            os.path.join(results_dir, f"{problem}_per_niche_success.csv"), index=False
        )

        niche_success = compile_niche_success(
            problem_df, num_parents=2, per_niche=False
        )

        success_cols = [
            "offspring_added_rate",
            "n_offspring_rank_mean",
            "offspring_dist_to_ideal_mean",
            "num_offspring",
        ]

        # save niche_success to csv
        niche_success.to_csv(
            os.path.join(results_dir, f"{problem}_niche_success.csv"), index=False
        )

        for success_col in success_cols:
            plot_niche_success(
                niche_success,
                max_niche=num_niches,
                output_dir=results_dir,
                file_prefix=problem,
                success_col=success_col,
            )

        distance_success = compile_distance_success(
            problem_df, num_parents=2, ideal_point=problem_dict["ideal_point"]
        )
        success_cols = [
            "offspring_dist_to_ideal",
            "n_offspring_rank",
            "ncrit_offspring_rank",
        ]

        # save distance_success to csv
        distance_success.to_csv(
            os.path.join(results_dir, f"{problem}_distance_success.csv"), index=False
        )
        pred_data = distance_success.copy()
        additional_cols = [
            "offspring_niche",
            "parent_0_niche",
            "parent_1_niche",
            "archive_stats_before_coverage",
            "archive_stats_before_unbounded_hv",
            "archive_stats_before_max_rank",
            "archive_stats_before_mean_rank",
            "archive_stats_after_coverage",
            "archive_stats_after_unbounded_hv",
            "archive_stats_after_max_rank",
            "archive_stats_after_mean_rank",
        ]
        pred_data = pred_data.merge(
            problem_df[additional_cols], how="left", left_index=True, right_index=True
        )
        pred_data.to_csv(
            os.path.join(results_dir, f"{problem}_pred_data.csv"), index=False
        )

        for success_col in success_cols:
            for dist_col in ["parent_x_distance", "parent_y_distance", "parent_angle"]:
                plot_distance_success_1d(
                    distance_success,
                    output_dir=results_dir,
                    dist_col=dist_col,
                    success_col=success_col,
                    file_prefix=problem,
                )
            plot_distance_success_2d(
                distance_success,
                output_dir=results_dir,
                x_axis="parent_x_distance",
                y_axis="parent_y_distance",
                success_col=success_col,
                file_prefix=problem,
            )
            plot_distance_success_2d(
                distance_success,
                output_dir=results_dir,
                x_axis="parent_x_distance",
                y_axis="parent_angle",
                success_col=success_col,
                file_prefix=problem,
            )


if __name__ == "__main__":
    results_dict = build_results_dict()
    # save results_dict to json file
    with open("results_dict.json", "w") as f:
        json.dump(results_dict, f, indent=4, default=str)
    postprocess_results(results_dict)
