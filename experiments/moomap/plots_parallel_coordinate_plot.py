import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import parallel_coordinates
from pathlib import Path


def plot_parallel_coordinate_plot(
    df: pd.DataFrame,
    class_column: str,
    features: list[str],
    title: str = "Parallel Coordinate Plot",
    save_path: str | Path | None = None,
):
    """
    Plots a parallel coordinate plot for the given dataframe.

    Args:
        df (pd.DataFrame): The dataframe containing the data to plot.
        features (list[str]): The list of feature names to include in the plot.
        target (str): The name of the target variable to color the lines by.
        title (str, optional): The title of the plot. Defaults to "Parallel Coordinate Plot".
        save_path (str | Path | None, optional): The path to save the plot. If None, the plot is shown. Defaults to None.
    """
    plt.figure(figsize=(12, 6))
    parallel_coordinates(df, class_column=class_column, cols=features)
    plt.title(title)
    plt.xlabel("Features")
    plt.ylabel("Values")
    plt.xticks(rotation=45)
    plt.grid()
    if save_path:
        plt.savefig(save_path, bbox_inches="tight")
        print(f"Plot saved to {save_path}")
