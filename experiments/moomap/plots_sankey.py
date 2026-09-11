from matplotlib.patches import PathPatch
from matplotlib.path import Path
import matplotlib.pyplot as plt
import pandas as pd
import os


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
) -> None:
    """
    Sankey flow plot for the top source and sink flows.
    """

    # preprocess df
    df = source_df.copy()

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
    df_grouped = df_selected.copy()

    # ---------------------------------------------------------
    # Nodes
    # ---------------------------------------------------------
    source_cols = df_grouped[source_col].unique().tolist()
    sink_cols = sorted(df_grouped[sink_col].unique())

    # Assign vertical positions
    # Top = 1, bottom = 0
    source_y = {
        source: 1 - (i + 1) / (len(source_cols) + 1)
        for i, source in enumerate(source_cols)
    }

    sink_y = {
        sink: 1 - (i + 1) / (len(sink_cols) + 1) for i, sink in enumerate(sink_cols)
    }

    # ---------------------------------------------------------
    # Scale flow widths
    # ---------------------------------------------------------
    max_flow = df_grouped[value_col].max()

    # Maximum vertical thickness of a flow

    # ---------------------------------------------------------
    # Figure
    # ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(14, 10))

    # ---------------------------------------------------------
    # Draw flows
    # ---------------------------------------------------------
    for _, row in df_grouped.iterrows():
        source = row[source_col]
        sink = row[sink_col]
        value = row[value_col]

        y0 = source_y[source]
        y1 = sink_y[sink]

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
    min_flow = df_grouped[value_col].min()
    max_flow = df_grouped[value_col].max()
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

    for source, y in source_y.items():
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
            source,
            ha="right",
            va="center",
            fontsize=9,
        )

    for sink, y in sink_y.items():
        ax.add_patch(
            plt.Rectangle(
                (right_x - node_width / 2, y - node_height / 2),
                node_width,
                node_height,
                color=sink_colour,
            )
        )

        ax.text(
            right_x + 0.025,
            y,
            str(sink),
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
        source_col,
        ha="center",
        va="bottom",
        fontsize=13,
        fontweight="bold",
    )

    ax.text(
        right_x,
        1.03,
        sink_col,
        ha="center",
        va="bottom",
        fontsize=13,
        fontweight="bold",
    )

    ax.set_title(
        f"{file_prefix}: Top {num_top_source_flows} {source_col} and Top {num_top_sink_flows} {sink_col} Flows",
        fontsize=15,
    )

    ax.axis("off")

    plt.tight_layout()

    file_name = f"{file_prefix}_sankey_{source_col}to{sink_col}_{value_col}.pdf"
    fig.savefig(
        os.path.join(
            output_dir,
            file_name,
        ),
        bbox_inches="tight",
    )

    plt.close(fig)
