import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_scatter(
    source_df: pd.DataFrame,
    output_dir: str,
    x_col: str,
    y_col: str,
    hue_col: str,
    file_prefix: str = "",
    cmap_str: str = "viridis",
) -> None:
    """
    Plot the success rate of offspring based on the distance between parents.
    """
    # filter rows with nan values in success_col

    df = source_df.copy()
    df = df.dropna(subset=[hue_col])

    norm = mpl.colors.Normalize(vmin=min(df[hue_col]), vmax=max(df[hue_col]))
    cmap = plt.get_cmap(cmap_str)

    plt.figure(figsize=(10, 8))
    ax = sns.scatterplot(
        data=df,
        x=x_col,
        y=y_col,
        hue=hue_col,
        palette=cmap,
        hue_norm=norm,
        alpha=0.7,
        legend=False,
    )
    # Create continuous colorbar
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    plt.colorbar(sm, ax=ax, label=hue_col)

    title = f"{file_prefix}: {hue_col} by {x_col} vs {y_col}"
    if len(source_df) > len(df):
        # We dropped some rows due to NaN values in the hue_col, so we should indicate that in the title
        title += " (dropped NaN)"
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    # plt.colorbar(label=success_col)
    plt.tight_layout()
    plt.savefig(
        os.path.join(
            output_dir, f"{file_prefix}_scatter_{x_col}vs{y_col}_{hue_col}.pdf"
        )
    )
    plt.close()
