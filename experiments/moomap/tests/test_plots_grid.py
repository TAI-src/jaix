import pandas as pd
import numpy as np
from plots_grid import plot_grid
import os


def test_plot_grid(tmp_path):
    # Create a sample DataFrame for testing
    df = pd.DataFrame(
        {
            "x": np.random.randint(0, 100, size=70),
            "y": np.random.randint(0, 100, size=70),
            "z": np.random.rand(70),
        }
    )
    df_grouped = (
        df.groupby(["x", "y"])
        .agg(z_mean=("z", "mean"), z_count=("z", "count"))
        .reset_index()
    )

    for col in ["z_mean", "z_count"]:
        assert col in df_grouped.columns

        file_path = plot_grid(
            source_df=df_grouped,
            max_grid=100,
            grid_colx="x",
            grid_coly="y",
            hue_col=col,
            output_dir=tmp_path,
            file_prefix="test_",
        )
        assert file_path.endswith(".pdf")
        assert os.path.exists(file_path)
        assert "xvsy" in file_path
        assert col in file_path
