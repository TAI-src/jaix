import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path


def select_top_flows(
    df: pd.DataFrame,
    source_col: str,
    sink_col: str,
    value_col: str,
    num_top_source_flows: int = 20,
    num_top_sink_flows: int = 5,
) -> pd.DataFrame:
    """
    Select the top source and sink flows from the DataFrame.
    """

    # Aggregate flows
    df_grouped = df.groupby([source_col, sink_col])[value_col].sum().reset_index()
    # Total value per source_col to determine the most important sources
    source_totals = (
        df_grouped.groupby(source_col)[value_col].sum().sort_values(ascending=False)
    )

    # Select the most important sources
    top_sources = source_totals.head(num_top_source_flows).index
    df_selected = df_grouped[df_grouped[source_col].isin(top_sources)].copy()

    # For each top source, keep its top K sink flows
    df_selected = (
        df_selected.sort_values(
            [source_col, value_col],
            ascending=[True, False],
        )
        .groupby(source_col)
        .head(num_top_sink_flows)
    )
    return df_selected


def get_node_order(
    df: pd.DataFrame, source_col: str, sink_col: str
) -> tuple[list[str], list[str]]:
    """
    Get the order of nodes for the Sankey plot.
    """

    source_cols = df[source_col].unique().tolist()
    sink_cols = sorted(df[sink_col].unique())
    return source_cols, sink_cols


def append_node_positions(
    df: pd.DataFrame,
    source_col: str,
    sink_col: str,
) -> pd.DataFrame:
    """
    Append node positions to the DataFrame for the Sankey plot.
    """

    source_cols, sink_cols = get_node_order(df, source_col, sink_col)

    # Assign vertical positions
    # Top = 1, bottom = 0
    source_y = {
        source: 1 - (i + 1) / (len(source_cols) + 1)
        for i, source in enumerate(source_cols)
    }

    sink_y = {
        sink: 1 - (i + 1) / (len(sink_cols) + 1) for i, sink in enumerate(sink_cols)
    }

    df["source_y"] = df[source_col].map(source_y)
    df["sink_y"] = df[sink_col].map(sink_y)

    return df


def append_width_scaling(
    df: pd.DataFrame,
    value_col: str,
    max_width: float = 0.025,
    min_width: float = 0.003,
) -> pd.DataFrame:
    """
    Append width scaling to the DataFrame for the Sankey plot.
    """

    max_flow = df[value_col].max()
    df["width"] = df[value_col].apply(
        lambda x: max(min_width, max_width * x / max_flow)
    )
    return df


def plot_specifications(
    df: pd.DataFrame,
    source_col: str,
    sink_col: str,
    value_col: str,
    max_width: float = 0.025,
    min_width: float = 0.003,
) -> pd.DataFrame:
    """
    Prepare the DataFrame for plotting by appending node positions.
    """

    df = append_node_positions(
        df,
        source_col=source_col,
        sink_col=sink_col,
    )
    df = append_width_scaling(
        df,
        value_col=value_col,
        max_width=max_width,
        min_width=min_width,
    )

    return df


def get_control_points(source_x, sink_x, bezier_control_factor):
    """
    Compute the control points for a cubic Bezier curve between source and sink.
    """

    dx = sink_x - source_x

    control_x1 = source_x + bezier_control_factor * dx
    control_x2 = sink_x - bezier_control_factor * dx

    return control_x1, control_x2


def compute_bezier_path(
    source_x: float,
    source_y: float,
    sink_x: float,
    sink_y: float,
    width: float,
    control_x1: float,
    control_x2: float,
) -> tuple[list[tuple[float, float]], list[np.uint8]]:
    """
    Compute the vertices and codes for a cubic Bezier path between source and sink.
    """
    # Upper boundary
    verts_top = [
        (source_x, source_y + width / 2),
        (control_x1, source_y + width / 2),
        (control_x2, sink_y + width / 2),
        (sink_x, sink_y + width / 2),
    ]

    # Lower boundary
    verts_bottom = [
        (sink_x, sink_y - width / 2),
        (control_x2, sink_y - width / 2),
        (control_x1, source_y - width / 2),
        (source_x, source_y - width / 2),
    ]

    vertices = verts_top + verts_bottom + [(source_x, source_y + width / 2)]

    codes = (
        [Path.MOVETO]
        + [Path.CURVE4] * 3
        + [Path.LINETO]
        + [Path.CURVE4] * 3
        + [Path.CLOSEPOLY]
    )

    return vertices, codes


def add_legend(
    ax,
    df,
    value_col,
    flow_colour,
    legend_x,
    legend_y,
    legend_spacing,
    max_width,
    min_width,
):
    """
    Add a legend to the Sankey plot indicating the flow values.
    """
    min_flow = df[value_col].min()
    max_flow = df[value_col].max()
    mid_flow = (min_flow + max_flow) / 2

    legend_values = [min_flow, mid_flow, max_flow]

    ax.text(
        legend_x,
        legend_y + 0.10,
        value_col,
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
            alpha=0.5,
            solid_capstyle="butt",
        )

        ax.text(
            legend_x + 0.17,
            y,
            f"{value:,.0f}",
            va="center",
            fontsize=9,
        )


def add_node(
    ax,
    x,
    y,
    node_width,
    node_height,
    color,
    label,
    label_offset=0.025,
    label_align="left",
):
    """
    Add a node to the Sankey plot.
    """
    ax.add_patch(
        Rectangle(
            (x - node_width / 2, y - node_height / 2),
            node_width,
            node_height,
            color=color,
        )
    )
    ax.text(
        x + label_offset if label_align == "left" else x - label_offset,
        y,
        label,
        ha=label_align,
        va="center",
        fontsize=9,
    )


def plot_sankey_flows(
    source_df: pd.DataFrame,
    output_dir: str,
    source_col: str,
    sink_col: str,
    value_col: str,
    file_prefix: str = "",
    num_top_source_flows: int = 20,
    num_top_sink_flows: int = 5,
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
    sink_colour: str = "darkorange",
    legend_x=0.35,
    legend_y=0.04,
    legend_spacing=0.035,
) -> str:
    """
    Sankey flow plot for the top source and sink flows.
    """

    # preprocess df
    df = select_top_flows(
        source_df,
        source_col=source_col,
        sink_col=sink_col,
        value_col=value_col,
        num_top_source_flows=num_top_source_flows,
        num_top_sink_flows=num_top_sink_flows,
    )

    df = plot_specifications(
        df,
        source_col=source_col,
        sink_col=sink_col,
        value_col=value_col,
        max_width=max_width,
        min_width=min_width,
    )

    fig, ax = plt.subplots(figsize=(14, 10))
    control_x1, control_x2 = get_control_points(
        source_x=left_x,
        sink_x=right_x,
        bezier_control_factor=bezier_control_factor,
    )
    for _, row in df.iterrows():
        # Draw the flow as a Bezier curve
        vertices, codes = compute_bezier_path(
            source_x=left_x,
            source_y=row["source_y"],
            sink_x=right_x,
            sink_y=row["sink_y"],
            width=row["width"],
            control_x1=control_x1,
            control_x2=control_x2,
        )

        patch = PathPatch(
            Path(vertices, codes),
            facecolor=flow_colour,
            edgecolor="none",
            alpha=alpha,
        )

        ax.add_patch(patch)

    # Go through distinct source and sink nodes to add them to the plot
    source_y = df.groupby(source_col)["source_y"].first().to_dict()
    sink_y = df.groupby(sink_col)["sink_y"].first().to_dict()

    for source, y in source_y.items():
        # Add the source node
        add_node(
            ax,
            x=left_x,
            y=y,
            node_width=node_width,
            node_height=node_height,
            color=source_colour,
            label=source,
            label_offset=node_width,
            label_align="right",
        )

    for sink, y in sink_y.items():
        ## Add the sink node
        add_node(
            ax,
            x=right_x,
            y=y,
            node_width=node_width,
            node_height=node_height,
            color=sink_colour,
            label=sink,
            label_offset=node_width,
            label_align="left",
        )

    # Add legend
    add_legend(
        ax=ax,
        df=df,
        value_col=value_col,
        flow_colour=flow_colour,
        legend_x=legend_x,
        legend_y=legend_y,
        legend_spacing=legend_spacing,
        max_width=max_width,
        min_width=min_width,
    )
    # Source label
    ax.text(
        left_x,
        1.03,
        source_col,
        ha="center",
        va="bottom",
        fontsize=13,
        fontweight="bold",
    )
    # Sink label
    ax.text(
        right_x,
        1.03,
        sink_col,
        ha="center",
        va="bottom",
        fontsize=13,
        fontweight="bold",
    )

    # ---------------------------------------------------------
    # Formatting
    # ---------------------------------------------------------
    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(0, 1)

    ax.set_title(
        f"{file_prefix}: Top {num_top_source_flows} {source_col} and Top {num_top_sink_flows} {sink_col} Flows",
        fontsize=15,
    )

    ax.axis("off")

    plt.tight_layout()

    file_name = f"{file_prefix}_sankey_{source_col}to{sink_col}_{value_col}.pdf"
    file_path = os.path.join(output_dir, file_name)
    fig.savefig(
        file_path,
        bbox_inches="tight",
    )

    plt.close(fig)
    return file_path
