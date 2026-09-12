from plots_scatter import plot_scatter
import pandas as pd
import numpy as np
import os


def test_plot_scatter(tmp_path):
    df = pd.DataFrame(
        {
            "x": np.random.normal(10, 2, 100),
            "y": np.random.normal(20, 5, 100),
            "z": np.random.uniform(0, 100, 100),
        }
    )

    plot_file = plot_scatter(
        df,
        tmp_path,
        x_col="x",
        y_col="y",
        hue_col="z",
        file_prefix="test_",
    )
    assert plot_file.endswith(".pdf")
    assert os.path.exists(plot_file)
    assert "xvsy" in plot_file
    assert "z" in plot_file
