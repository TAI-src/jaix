import numpy as np
import pandas as pd
import os
import pytest
from mating_kernel.plots.heatmap import plot_heatmap

test_data_regular = [
    {
        "setting": "setting",
        "pid": "pid",
        "cv_score_mean": 0.5,
        "cv_score_std": 0.1,
    },
    {
        "setting": "setting",
        "pid": "pid2",
        "cv_score_mean": 0.6,
        "cv_score_std": 0.2,
    },
    {
        "setting": "setting2",
        "pid": "pid",
        "cv_score_mean": 0.7,
        "cv_score_std": 0.3,
    },
    {
        "setting": "setting2",
        "pid": "pid2",
        "cv_score_mean": 0.8,
        "cv_score_std": 0.4,
    },
]

test_data_outliers = [
    {
        "setting": "setting",
        "pid": "pid",
        "cv_score_mean": 0.5,
        "cv_score_std": 0.1,
    },
    {
        "setting": "setting",
        "pid": "pid2",
        "cv_score_mean": np.nan,
        "cv_score_std": 285724.3218554493,
    },
]


@pytest.mark.parametrize("test_data", [test_data_regular, test_data_outliers])
def test_plot_heatmap(tmp_path, test_data):

    df = pd.DataFrame(test_data)
    file_path = plot_heatmap(
        source_df=df,
        max_grid=None,
        grid_colx="pid",
        grid_coly="setting",
        hue_col="cv_score_mean",
        output_dir=tmp_path,
        file_prefix=str(np.random.randint(0, 10000)),
        cmap_pal="viridis",
        annot=True,
    )
    assert file_path.endswith(".pdf")
    assert os.path.exists(file_path)
