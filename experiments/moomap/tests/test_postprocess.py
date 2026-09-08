import os
import pandas as pd
import numpy as np


from postprocess import (
    get_success_cols,
    compile_niche_success,
    compile_distance_success,
    plot_niche_success,
    plot_distance_success_2d,
    plot_distance_success_1d,
)


test_data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_results.csv")


def get_test_data():
    df = pd.read_csv(test_data)
    gen_dict = {}
    for gen in range(10):
        gen_df = df[df["generation"] == gen]
        max_rank = gen_df["offspring_rank"].max()
        added_rows = gen_df[gen_df["offspring_added"] == True]
        crit_rank = added_rows["offspring_rank"].max()
        gen_dict[gen] = {
            "max_rank": max_rank,
            "crit_rank": crit_rank,
            "n_added": len(added_rows),
            "add_rate": len(added_rows) / len(gen_df),
            "min_rank": gen_df["offspring_rank"].min(),
            "mean_rank": gen_df["offspring_rank"].mean(),
        }
    return df, gen_dict


def test_get_success_cols():
    df, gen_dict = get_test_data()
    success_cols = get_success_cols(df)
    assert set(success_cols) == {
        "offspring_added",
        "offspring_rank",
        "n_offspring_rank",
        "ncrit_offspring_rank",
    }
    assert "n_offspring_rank" in df.columns
    assert "ncrit_offspring_rank" in df.columns
    assert "critical_rank" in df.columns
    for gen in range(10):
        row_idx = gen * 100
        row_data = df.loc[row_idx]
        assert row_data["generation"] == gen
        assert row_data["critical_rank"] == gen_dict[gen]["crit_rank"]

    # check n_offspring_rank is between 0 and 1
    assert df["n_offspring_rank"].between(0, 1).all()
    # check that where offspring_rank = 0,
    assert df[df["offspring_rank"] == 0]["n_offspring_rank"].eq(0).all()
    assert df[df["ncrit_offspring_rank"] == 0]["ncrit_offspring_rank"].eq(0).all()
    # check that where offspring_rank = max_rank, n_offspring_rank = 1
    # check that where offspring_rank = crit_rank, ncrit_offspring_rank = 1
    for gen in range(10):
        max_rank = gen_dict[gen]["max_rank"]
        assert (
            df[(df["generation"] == gen) & (df["offspring_rank"] == max_rank)][
                "n_offspring_rank"
            ]
            .eq(1)
            .all()
        )
        crit_rank = gen_dict[gen]["crit_rank"]
        if crit_rank > 0:
            assert (
                df[(df["generation"] == gen) & (df["offspring_rank"] == crit_rank)][
                    "ncrit_offspring_rank"
                ]
                .eq(1)
                .all()
            )
            if max_rank > crit_rank:
                # check that where offspring_rank > crit_rank, ncrit_offspring_rank > 1
                assert (
                    df[(df["generation"] == gen) & (df["offspring_rank"] > crit_rank)][
                        "ncrit_offspring_rank"
                    ]
                    .gt(1)
                    .all()
                )


def test_compile_niche_success():
    df, gen_dict = get_test_data()
    df1 = df[(df["parent_0_niche"] == 6) & (df["parent_1_niche"] == 34)][
        "offspring_added"
    ]
    df2 = df[(df["parent_0_niche"] == 34) & (df["parent_1_niche"] == 6)][
        "offspring_added"
    ]
    expected_num_offspring = len(df1) + len(df2)
    expected_added_rate = (df1.sum() + df2.sum()) / expected_num_offspring
    niche_success = compile_niche_success(df)
    # check that the number of rows is smaller than the number of unique parent niche pairs
    # This is because the order of the parent niches does not matter, so (niche_a, niche_b) is the same as (niche_b, niche_a)
    # however, the original df should not be changed
    assert len(niche_success) < len(df.groupby(["parent_0_niche", "parent_1_niche"]))
    # check that the offspring_added_rate is between 0 and 1
    assert niche_success["offspring_added_rate"].between(0, 1).all()

    niche_success_01 = niche_success[
        (niche_success["parent_0_niche"] == 6) & (niche_success["parent_1_niche"] == 34)
    ]
    assert len(niche_success_01) == 1
    assert niche_success_01["num_offspring"].iloc[0] == expected_num_offspring
    assert niche_success_01["offspring_added_rate"].iloc[0] == expected_added_rate


def test_plot_niche_success(tmp_path):
    df, _ = get_test_data()
    niche_success = compile_niche_success(df)
    os.makedirs(tmp_path, exist_ok=True)
    plot_niche_success(
        niche_success,
        100,
        tmp_path,
        success_col="offspring_added_rate",
        file_prefix="test_",
    )


def test_compile_distance_success():

    df, gen_dict = get_test_data()
    parent_0_x = df["parent_0_x"].iloc[0]
    parent_0_y = df["parent_0_y"].iloc[0]
    parent_1_x = df["parent_1_x"].iloc[0]
    parent_1_y = df["parent_1_y"].iloc[0]
    # convert string representation of arrays to actual arrays
    p0x = np.fromstring(parent_0_x.strip("[]"), sep=" ")
    p0y = np.fromstring(parent_0_y.strip("[]"), sep=" ")
    p1x = np.fromstring(parent_1_x.strip("[]"), sep=" ")
    p1y = np.fromstring(parent_1_y.strip("[]"), sep=" ")
    expected_x_distance = np.linalg.norm(p0x - p1x)
    expected_y_distance = np.linalg.norm(p0y - p1y)
    distance_success = compile_distance_success(df)
    # check that the number of rows is equal to the number of rows in the original df
    assert len(distance_success) == len(df)
    # check that the offspring_added_rate is between 0 and 1
    assert distance_success["offspring_added"].between(0, 1).all()
    # check that the parent distances are correct based on the first one
    assert np.isclose(
        distance_success["parent_x_distance"].iloc[0], expected_x_distance
    )
    assert np.isclose(
        distance_success["parent_y_distance"].iloc[0], expected_y_distance
    )
    assert "parent_0_dist_to_ideal" in distance_success.columns
    assert "parent_1_dist_to_ideal" in distance_success.columns
    assert "offspring_dist_to_ideal" in distance_success.columns
    assert "n_offspring_rank" in distance_success.columns
    assert "ncrit_offspring_rank" in distance_success.columns


def test_plot_distance_success(tmp_path):
    out_dir = tmp_path
    os.makedirs(out_dir, exist_ok=True)
    df, _ = get_test_data()
    distance_success = compile_distance_success(df)
    num_rows = len(distance_success)
    plot_distance_success_2d(
        distance_success,
        out_dir,
        success_col="offspring_dist_to_ideal",
        file_prefix="test_",
    )
    assert len(distance_success) == num_rows
    plot_distance_success_1d(
        distance_success,
        out_dir,
        success_col="offspring_dist_to_ideal",
        file_prefix="test_",
    )
