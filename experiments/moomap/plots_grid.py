import os

import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns


def plot_grid(
    source_df: pd.DataFrame,
    max_grid: int,
    grid_colx: str,
    grid_coly: str,
    success_col: str,
    output_dir: str,
    file_prefix: str = "",
    cmap_pal: str = "viridis",
) -> None:
    """
    Plot the success rate of each niche pair in the data.
    """
    grid_data = source_df.copy()
    # create a pivot table with the mean success rate for each grid pair
    # grid labels should be from 0 to max_grid, so we can create a pivot table for the complete grid
    idx = pd.MultiIndex.from_product(
        [range(max_grid), range(max_grid)], names=[grid_colx, grid_coly]
    )
    grid_data = grid_data.set_index([grid_colx, grid_coly]).reindex(idx).reset_index()
    pivot_table = grid_data.pivot(
        index=grid_colx, columns=grid_coly, values=success_col
    )
    cmap = sns.color_palette(cmap_pal, as_cmap=True)
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
    plt.title(f"{file_prefix}: {success_col} by {grid_colx} vs {grid_coly}")
    plt.xlabel(grid_colx)
    plt.ylabel(grid_coly)
    plt.tight_layout()
    file_path = os.path.join(
        output_dir, f"{file_prefix}_grid_{grid_colx}vs{grid_coly}_{success_col}.pdf"
    )
    plt.savefig(file_path)
    plt.close()
