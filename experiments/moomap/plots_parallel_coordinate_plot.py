import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D


def plot_pcp(
    df: pd.DataFrame,
    line_col_name: str,
    class_col_name: str | None = None,
    linestyles: dict[str, str] | None = None,
    features: list[str] | None = None,
    cmap_name: str = "tab10",
    save_path: str | None = None,
):

    fig, ax = plt.subplots(figsize=(12, 6))

    # One colour per unique value in line_col
    names = df[line_col_name].unique()
    cmap = plt.colormaps[cmap_name]

    colors = [cmap(i / max(len(names) - 1, 1)) for i in range(len(names))]

    name_colors = dict(zip(names, colors))

    if features is None:
        features = list(
            df.columns.drop(
                [line_col_name, class_col_name] if class_col_name else [line_col_name]
            )
        )
    x = range(len(features))

    for _, row in df.iterrows():
        if linestyles is None:
            ax.plot(
                x,
                row[features].values,
                color=name_colors[row[line_col_name]],
                alpha=0.8,
            )
        else:
            assert (
                class_col_name is not None
            ), "class_col_name must be provided if linestyles is provided"

            ax.plot(
                x,
                row[features].values,
                color=name_colors[row[line_col_name]],
                linestyle=linestyles.get(row[class_col_name], "-"),
                alpha=0.8,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(features, rotation=45, ha="right")
    ax.grid(
        axis="both",
        linestyle=":",
        linewidth=0.8,
        alpha=0.5,
    )

    # name legend only
    handles = [
        Line2D(
            [0],
            [0],
            color=name_colors[name],
            lw=2,
            label=name,
        )
        for name in names
    ]

    ax.legend(handles=handles, title=line_col_name)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
