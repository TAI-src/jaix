import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_grid(
    source_df: pd.DataFrame,
    max_grid: int | None,
    grid_colx: str,
    grid_coly: str,
    hue_col: str,
    output_dir: str,
    file_prefix: str = "",
    cmap_pal: str = "viridis",
) -> str:
    """
    Plots a grid heatmap using hue_col as the color intensity for each grid cell defined by grid_colx and grid_coly.
    """
    grid_data = source_df.copy()
    # create a pivot table using grid_colx and grid_coly as indices and hue_col as values
    # grid labels should be from 0 to max_grid, so we can create a pivot table for the complete grid
    if max_grid is not None:
        idx = pd.MultiIndex.from_product(
            [range(max_grid), range(max_grid)], names=[grid_colx, grid_coly]
        )
        grid_data = (
            grid_data.set_index([grid_colx, grid_coly]).reindex(idx).reset_index()
        )
    pivot_table = grid_data.pivot(index=grid_colx, columns=grid_coly, values=hue_col)
    cmap = sns.color_palette(cmap_pal, as_cmap=True)
    cmap.set_bad("lightgray")
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        pivot_table,
        annot=False,
        fmt=".2f",
        cmap=cmap,
        cbar_kws={"label": hue_col},
    )
    plt.gca().invert_yaxis()
    plt.title(f"{file_prefix}: {hue_col} by {grid_colx} vs {grid_coly}")
    plt.xlabel(grid_coly)
    plt.ylabel(grid_colx)
    plt.tight_layout()
    file_name = f"{file_prefix}_grid_{grid_colx}vs{grid_coly}_{hue_col}.pdf"
    file_path = os.path.join(output_dir, file_name)
    plt.savefig(file_path)
    plt.close()
    return file_path
