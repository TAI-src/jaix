import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from postprocess import (
    compile_distance_success,
    compile_niche_success,
    compile_pred_data,
    get_success_cols,
    plot_sankey,
    postprocess_results,
    run_postprocess,
)
from utils_read import get_nsga3x_results

data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

test_data = os.path.join(data_path, "test_results.csv")
ENABLE_PLOT_TESTS = False


def get_test_data():
    df = pd.read_csv(test_data)
    gen_dict = {}
    for gen in range(10):
        gen_df = df[df["generation"] == gen]
        max_rank = gen_df["offspring_rank"].max()
        added_rows = gen_df[gen_df["offspring_added"] == True]
        removed_rows = gen_df[gen_df["offspring_added"] == False]
        crit_rank = removed_rows["offspring_rank"].min()
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
        "offspring_dist_to_ideal",
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


@pytest.mark.parametrize("per_niche", [True, False])
def test_compile_niche_success(tmp_path, per_niche):
    df, _gen_dict = get_test_data()
    df1 = df[(df["parent_0_niche"] == 6) & (df["parent_1_niche"] == 34)][
        "offspring_added"
    ]
    df2 = df[(df["parent_0_niche"] == 34) & (df["parent_1_niche"] == 6)][
        "offspring_added"
    ]
    expected_num_offspring = len(df1) + len(df2)
    expected_added_rate = (df1.sum() + df2.sum()) / expected_num_offspring
    niche_success, niche_file = compile_niche_success(
        df, per_niche=per_niche, output_dir=tmp_path, file_prefix="test_niche_success"
    )
    assert niche_file is not None
    assert os.path.exists(niche_file)
    # check that the number of rows is smaller than the number of unique parent niche pairs
    # This is because the order of the parent niches does not matter, so (niche_a, niche_b) is the same as (niche_b, niche_a)
    # however, the original df should not be changed
    if per_niche:
        assert len(niche_success) < len(
            df.groupby(["parent_0_niche", "parent_1_niche", "offspring_niche"])
        )
    else:
        assert len(niche_success) < len(
            df.groupby(["parent_0_niche", "parent_1_niche"])
        )
    # check that the offspring_added_rate is between 0 and 1
    assert niche_success["offspring_added_rate"].between(0, 1).all()

    niche_success_01 = niche_success[
        (niche_success["parent_0_niche"] == 6) & (niche_success["parent_1_niche"] == 34)
    ]
    if per_niche:
        # Two offspring niches for this niche pair
        assert len(niche_success_01) == 2
        assert sum(niche_success_01["num_offspring"]) == expected_num_offspring
        assert np.isclose(
            sum(
                niche_success_01["offspring_added_rate"]
                * niche_success_01["num_offspring"]
            )
            / expected_num_offspring,
            expected_added_rate,
        )
    else:
        assert len(niche_success_01) == 1
        assert niche_success_01["num_offspring"].iloc[0] == expected_num_offspring
        assert niche_success_01["offspring_added_rate"].iloc[0] == expected_added_rate


def test_compile_distance_success(tmp_path):

    df, _gen_dict = get_test_data()
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
    distance_success, dist_file = compile_distance_success(
        df,
        ideal_point=np.array([0, 0]),
        output_dir=tmp_path,
        file_prefix="test_distance_success",
    )
    assert dist_file is not None
    assert os.path.exists(dist_file)
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
    assert "parent_angle" in distance_success.columns


def test_plot_sankey(tmp_path):
    df, _gen_dict = get_test_data()
    niche_success, _ = compile_niche_success(df, per_niche=True)
    # check that the function runs without error
    file_path = plot_sankey(
        niche_success, output_dir=tmp_path, file_prefix="test_sankey"
    )
    assert os.path.exists(file_path)


def test_compile_pred_data(tmp_path):
    df, gen_dict = get_test_data()
    distance_success, _ = compile_distance_success(df)
    pred_data, pred_file = compile_pred_data(
        distance_success, df, results_dir=tmp_path, problem="test_problem"
    )
    assert pred_file is not None
    assert os.path.exists(pred_file)
    # check that the number of rows is equal to the number of rows in the original df
    assert len(pred_data) == len(distance_success)
    # check that the columns are correct
    expected_cols = [
        "offspring_niche",
        "parent_0_niche",
        "parent_1_niche",
        "archive_stats_before_coverage",
        "archive_stats_before_unbounded_hv",
        "archive_stats_before_max_rank",
        "archive_stats_before_mean_rank",
        "archive_stats_after_coverage",
        "archive_stats_after_unbounded_hv",
        "archive_stats_after_max_rank",
        "archive_stats_after_mean_rank",
    ]
    for col in expected_cols:
        assert col in pred_data.columns


@pytest.mark.parametrize("skip_plot", [True, False])
def test_postprocess_results(tmp_path, skip_plot):
    if not ENABLE_PLOT_TESTS and not skip_plot:
        pytest.skip("Skipping plot tests because ENABLE_PLOT_TESTS is False")
    # Test that the postprocess_results function runs without error
    test_folder = Path(__file__).parent / "data" / "nsga3x_results"
    problem_ids = [0, 1]
    results_dict = get_nsga3x_results(test_folder, problem_ids=problem_ids)
    res = postprocess_results(results_dict, out_dir=tmp_path, skip_plots=skip_plot)
    assert isinstance(res, dict)
    assert res.keys() == set(problem_ids)
    for problem_id in problem_ids:
        assert problem_id in res
        assert len(res[problem_id]["data"]) == 4
        for data_file in res[problem_id]["data"].values():
            assert os.path.exists(data_file)
        if skip_plot:
            assert "plots" not in res[problem_id]
        else:
            assert "plots" in res[problem_id]
            assert len(res[problem_id]["plots"]) == 4
            sankey_file = res[problem_id]["plots"].pop("sankey")
            assert os.path.exists(sankey_file)
            for plot_type, plot_list in res[problem_id]["plots"].items():
                assert len(plot_list) >= 3  # minimum success cols is 3
                for plot_file in plot_list:
                    assert os.path.exists(plot_file)


def test_run_postprocess(tmp_path):
    # Test that the run_postprocess function runs without error
    test_folder = Path(__file__).parent / "data" / "nsga3x_results"
    problem_ids = [0, 1]
    res = run_postprocess(
        results_dir=str(test_folder),
        problem_ids=problem_ids,
        out_dir=tmp_path,
        skip_plots=True,
    )
    assert isinstance(res, dict)
    assert all(p in res for p in problem_ids)
    assert "meta" in res
    res_file = res["meta"]["res_file"]
    assert os.path.exists(res_file)
    for pid in problem_ids:
        value = res[pid]
        assert isinstance(value, dict)
        assert "data" in value
        assert "ideal_point" in value
