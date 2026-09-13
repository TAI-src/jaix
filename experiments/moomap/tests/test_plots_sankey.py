import os

import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.patches import Path

from plots_sankey import (
    add_legend,
    add_node,
    append_node_positions,
    append_width_scaling,
    compute_bezier_path,
    get_control_points,
    get_node_order,
    plot_sankey_flows,
    plot_specifications,
    select_top_flows,
)


def test_plot_sankey_flows(tmp_path):

    # Create a sample DataFrame
    data = {
        "source": ["A", "A", "B", "B", "C", "C"],
        "sink": ["X", "Y", "X", "Z", "Y", "Z"],
        "value": [10, 5, 15, 20, 25, 30],
    }
    df = pd.DataFrame(data)

    # Call the function
    file_path = plot_sankey_flows(
        source_df=df,
        output_dir=tmp_path,
        source_col="source",
        sink_col="sink",
        value_col="value",
        file_prefix="test_",
        num_top_source_flows=2,
        num_top_sink_flows=1,
    )
    assert "test_" in file_path
    assert os.path.exists(file_path)
    assert file_path.endswith(".pdf")
    assert "sourcetosink" in file_path
    assert "value" in file_path


# Below tests are generated to cover the individual functions in plots_sankey.py. They are not exhaustive but provide a good starting point for testing the functionality of the module.


# ---------------------------------------------------------------------------
# Fixtures / test data
# ---------------------------------------------------------------------------


@pytest.fixture
def flow_df():
    """Small, deterministic flow dataset."""
    return pd.DataFrame(
        {
            "source": [
                "A",
                "A",
                "A",
                "B",
                "B",
                "C",
                "C",
            ],
            "sink": [
                "X",
                "Y",
                "Z",
                "X",
                "Y",
                "X",
                "Z",
            ],
            "value": [
                100,
                50,
                10,
                80,
                20,
                40,
                30,
            ],
        }
    )


# ---------------------------------------------------------------------------
# select_top_flows
# ---------------------------------------------------------------------------


def test_select_top_flows_aggregates_duplicate_flows():
    df = pd.DataFrame(
        {
            "source": ["A", "A", "A", "B"],
            "sink": ["X", "X", "Y", "X"],
            "value": [10, 20, 30, 40],
        }
    )

    result = select_top_flows(
        df,
        source_col="source",
        sink_col="sink",
        value_col="value",
    )

    # A -> X should have been aggregated from 10 + 20.
    row = result[(result["source"] == "A") & (result["sink"] == "X")]

    assert len(row) == 1
    assert row.iloc[0]["value"] == 30


def test_select_top_flows_selects_top_sources(flow_df):
    result = select_top_flows(
        flow_df,
        source_col="source",
        sink_col="sink",
        value_col="value",
        num_top_source_flows=2,
        num_top_sink_flows=5,
    )

    # Totals:
    # A = 160
    # B = 100
    # C = 70
    assert set(result["source"]) == {"A", "B"}


def test_select_top_flows_selects_top_sinks_per_source():
    df = pd.DataFrame(
        {
            "source": ["A", "A", "A", "B", "B"],
            "sink": ["X", "Y", "Z", "X", "Y"],
            "value": [100, 50, 10, 80, 20],
        }
    )

    result = select_top_flows(
        df,
        source_col="source",
        sink_col="sink",
        value_col="value",
        num_top_source_flows=2,
        num_top_sink_flows=1,
    )

    # Only the largest sink for each source should remain.
    assert set(zip(result["source"], result["sink"])) == {
        ("A", "X"),
        ("B", "X"),
    }


def test_select_top_flows_respects_both_limits(flow_df):
    result = select_top_flows(
        flow_df,
        source_col="source",
        sink_col="sink",
        value_col="value",
        num_top_source_flows=2,
        num_top_sink_flows=1,
    )

    assert result["source"].nunique() == 2
    assert len(result) == 2


# ---------------------------------------------------------------------------
# get_node_order
# ---------------------------------------------------------------------------


def test_get_node_order_preserves_source_order_and_sorts_sinks():
    df = pd.DataFrame(
        {
            "source": ["B", "A", "B", "C"],
            "sink": ["Z", "Y", "X", "Y"],
            "value": [1, 2, 3, 4],
        }
    )

    sources, sinks = get_node_order(
        df,
        source_col="source",
        sink_col="sink",
    )

    assert sources == ["B", "A", "C"]
    assert sinks == ["X", "Y", "Z"]


def test_get_node_order_returns_unique_nodes():
    df = pd.DataFrame(
        {
            "source": ["A", "A", "B", "B"],
            "sink": ["X", "X", "Y", "Y"],
        }
    )

    sources, sinks = get_node_order(
        df,
        source_col="source",
        sink_col="sink",
    )

    assert sources == ["A", "B"]
    assert sinks == ["X", "Y"]


# ---------------------------------------------------------------------------
# append_node_positions
# ---------------------------------------------------------------------------


def test_append_node_positions_adds_expected_columns(flow_df):
    result = append_node_positions(
        flow_df.copy(),
        source_col="source",
        sink_col="sink",
    )

    assert "source_y" in result.columns
    assert "sink_y" in result.columns


def test_append_node_positions_assigns_unique_source_positions(flow_df):
    result = append_node_positions(
        flow_df.copy(),
        source_col="source",
        sink_col="sink",
    )

    source_positions = result.groupby("source")["source_y"].first()

    assert source_positions.nunique() == 3

    # Positions should be between 0 and 1.
    assert source_positions.between(0, 1).all()


def test_append_node_positions_assigns_unique_sink_positions(flow_df):
    result = append_node_positions(
        flow_df.copy(),
        source_col="source",
        sink_col="sink",
    )

    sink_positions = result.groupby("sink")["sink_y"].first()

    assert sink_positions.nunique() == 3
    assert sink_positions.between(0, 1).all()


def test_append_node_positions_same_node_has_same_position(flow_df):
    result = append_node_positions(
        flow_df.copy(),
        source_col="source",
        sink_col="sink",
    )

    for _, group in result.groupby("source"):
        assert group["source_y"].nunique() == 1

    for _, group in result.groupby("sink"):
        assert group["sink_y"].nunique() == 1


# ---------------------------------------------------------------------------
# append_width_scaling
# ---------------------------------------------------------------------------


def test_append_width_scaling_adds_width_column():
    df = pd.DataFrame({"value": [10, 20, 30]})

    result = append_width_scaling(
        df,
        value_col="value",
    )

    assert "width" in result.columns


def test_append_width_scaling_max_value_gets_max_width():
    df = pd.DataFrame({"value": [10, 20, 30]})

    result = append_width_scaling(
        df,
        value_col="value",
        max_width=0.025,
        min_width=0.003,
    )

    assert result.loc[result["value"] == 30, "width"].iloc[0] == 0.025


def test_append_width_scaling_respects_min_width():
    df = pd.DataFrame({"value": [1, 100]})

    result = append_width_scaling(
        df,
        value_col="value",
        max_width=0.025,
        min_width=0.003,
    )

    assert result["width"].min() >= 0.003


def test_append_width_scaling_scales_proportionally():
    df = pd.DataFrame({"value": [25, 50, 100]})

    result = append_width_scaling(
        df,
        value_col="value",
        max_width=0.04,
        min_width=0.001,
    )

    assert result.loc[0, "width"] == 0.01
    assert result.loc[1, "width"] == 0.02
    assert result.loc[2, "width"] == 0.04


# ---------------------------------------------------------------------------
# plot_specifications
# ---------------------------------------------------------------------------


def test_plot_specifications_adds_plot_columns(flow_df):
    result = plot_specifications(
        flow_df.copy(),
        source_col="source",
        sink_col="sink",
        value_col="value",
    )

    assert "source_y" in result.columns
    assert "sink_y" in result.columns
    assert "width" in result.columns


def test_plot_specifications_preserves_original_data(flow_df):
    result = plot_specifications(
        flow_df.copy(),
        source_col="source",
        sink_col="sink",
        value_col="value",
    )

    pd.testing.assert_frame_equal(
        result[flow_df.columns],
        flow_df,
    )


# ---------------------------------------------------------------------------
# get_control_points
# ---------------------------------------------------------------------------


def test_get_control_points():
    control_x1, control_x2 = get_control_points(
        source_x=0.0,
        sink_x=1.0,
        bezier_control_factor=0.4,
    )

    assert control_x1 == 0.4
    assert control_x2 == 0.6


def test_get_control_points_with_factor_zero():
    control_x1, control_x2 = get_control_points(
        source_x=0.2,
        sink_x=0.8,
        bezier_control_factor=0,
    )

    assert control_x1 == 0.2
    assert control_x2 == 0.8


def test_get_control_points_with_factor_half():
    control_x1, control_x2 = get_control_points(
        source_x=0.0,
        sink_x=1.0,
        bezier_control_factor=0.5,
    )

    assert control_x1 == 0.5
    assert control_x2 == 0.5


# ---------------------------------------------------------------------------
# compute_bezier_path
# ---------------------------------------------------------------------------


def test_compute_bezier_path_returns_expected_number_of_vertices():
    vertices, codes = compute_bezier_path(
        source_x=0.0,
        source_y=0.5,
        sink_x=1.0,
        sink_y=0.6,
        width=0.1,
        control_x1=0.4,
        control_x2=0.6,
    )

    assert len(vertices) == 9
    assert len(codes) == 9


def test_compute_bezier_path_starts_and_ends_at_source():
    vertices, _codes = compute_bezier_path(
        source_x=0.0,
        source_y=0.5,
        sink_x=1.0,
        sink_y=0.6,
        width=0.1,
        control_x1=0.4,
        control_x2=0.6,
    )

    # First and final vertices should be the same point.
    assert vertices[0] == vertices[-1]


def test_compute_bezier_path_has_expected_path_codes():
    _vertices, codes = compute_bezier_path(
        source_x=0.0,
        source_y=0.5,
        sink_x=1.0,
        sink_y=0.6,
        width=0.1,
        control_x1=0.4,
        control_x2=0.6,
    )

    assert codes[0] == Path.MOVETO
    assert codes[1:4] == [Path.CURVE4] * 3
    assert codes[4] == Path.LINETO
    assert codes[5:8] == [Path.CURVE4] * 3
    assert codes[8] == Path.CLOSEPOLY


def test_compute_bezier_path_has_correct_width():
    width = 0.2

    vertices, _ = compute_bezier_path(
        source_x=0.0,
        source_y=0.5,
        sink_x=1.0,
        sink_y=0.5,
        width=width,
        control_x1=0.4,
        control_x2=0.6,
    )

    top_y = vertices[0][1]
    bottom_y = vertices[7][1]

    assert top_y == 0.5 + width / 2
    assert bottom_y == 0.5 - width / 2


# ---------------------------------------------------------------------------
# add_node
# ---------------------------------------------------------------------------


def test_add_node_adds_rectangle_and_label():
    fig, ax = plt.subplots()

    add_node(
        ax=ax,
        x=0.5,
        y=0.5,
        node_width=0.1,
        node_height=0.2,
        color="black",
        label="Test",
    )

    assert len(ax.patches) == 1
    assert len(ax.texts) == 1

    assert ax.texts[0].get_text() == "Test"

    plt.close(fig)


def test_add_node_uses_correct_dimensions():
    fig, ax = plt.subplots()

    add_node(
        ax=ax,
        x=0.5,
        y=0.5,
        node_width=0.1,
        node_height=0.2,
        color="red",
        label="Test",
    )

    rectangle = ax.patches[0]

    assert rectangle.get_width() == 0.1
    assert rectangle.get_height() == 0.2

    plt.close(fig)


def test_add_node_respects_label_alignment():
    fig, ax = plt.subplots()

    add_node(
        ax=ax,
        x=0.5,
        y=0.5,
        node_width=0.1,
        node_height=0.2,
        color="black",
        label="Test",
        label_offset=0.025,
        label_align="right",
    )

    text = ax.texts[0]

    assert text.get_ha() == "right"
    assert text.get_position() == (0.475, 0.5)

    plt.close(fig)


# ---------------------------------------------------------------------------
# add_legend
# ---------------------------------------------------------------------------


def test_add_legend_adds_title_values_and_lines():
    df = pd.DataFrame({"value": [10, 50, 100]})

    fig, ax = plt.subplots()

    add_legend(
        ax=ax,
        df=df,
        value_col="value",
        flow_colour="steelblue",
        legend_x=0.35,
        legend_y=0.04,
        legend_spacing=0.035,
        max_width=0.025,
        min_width=0.003,
    )

    # One title + three value labels.
    assert len(ax.texts) == 4

    # Three legend lines.
    assert len(ax.lines) == 3

    assert ax.texts[0].get_text() == "value"

    plt.close(fig)


def test_add_legend_uses_min_mid_and_max_values():
    df = pd.DataFrame({"value": [10, 50, 100]})

    fig, ax = plt.subplots()

    add_legend(
        ax=ax,
        df=df,
        value_col="value",
        flow_colour="steelblue",
        legend_x=0.35,
        legend_y=0.04,
        legend_spacing=0.035,
        max_width=0.025,
        min_width=0.003,
    )

    labels = [text.get_text() for text in ax.texts]

    assert "10" in labels
    assert "55" in labels
    assert "100" in labels

    plt.close(fig)


# ---------------------------------------------------------------------------
# plot_sankey_flows
# ---------------------------------------------------------------------------


def test_plot_sankey_flows_creates_pdf(tmp_path, flow_df):
    file_path = plot_sankey_flows(
        source_df=flow_df,
        output_dir=str(tmp_path),
        source_col="source",
        sink_col="sink",
        value_col="value",
        file_prefix="test",
    )

    assert file_path.endswith("test_sankey_sourcetosink_value.pdf")

    assert tmp_path.joinpath("test_sankey_sourcetosink_value.pdf").exists()


def test_plot_sankey_flows_returns_file_path(tmp_path, flow_df):
    result = plot_sankey_flows(
        source_df=flow_df,
        output_dir=str(tmp_path),
        source_col="source",
        sink_col="sink",
        value_col="value",
        file_prefix="example",
    )

    expected = tmp_path / "example_sankey_sourcetosink_value.pdf"

    assert result == str(expected)


def test_plot_sankey_flows_respects_top_flow_limits(tmp_path):
    df = pd.DataFrame(
        {
            "source": ["A"] * 4 + ["B"] * 4 + ["C"] * 4,
            "sink": ["W", "X", "Y", "Z"] * 3,
            "value": [
                100,
                90,
                80,
                70,
                60,
                50,
                40,
                30,
                20,
                10,
                5,
                1,
            ],
        }
    )

    result = plot_sankey_flows(
        source_df=df,
        output_dir=str(tmp_path),
        source_col="source",
        sink_col="sink",
        value_col="value",
        file_prefix="test",
        num_top_source_flows=2,
        num_top_sink_flows=2,
    )

    assert result.endswith("test_sankey_sourcetosink_value.pdf")
    assert tmp_path.joinpath("test_sankey_sourcetosink_value.pdf").exists()
