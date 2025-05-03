import wandb
import numpy as np
from collections import defaultdict
import json
import itertools
import pickle
import argparse
import pandas as pd


def bin_d(xi, xj):
    """
    Calculate the binary distance between two binary vectors.
    The binary distance is defined as the sum of absolute differences
    between the corresponding elements of the vectors.
    """
    return sum(abs(np.array(xi) - np.array(xj)))


def which2dec(x):
    """
    Convert a binary vector in "which" representation to its decimal representation.
    The binary vector is represented as a list of indices where the value is 1.
    """
    return [sum([2**i for i in xi]) for xi in x]


def which2bin(x, n):
    """
    Convert a binary vector in "which" representation to its binary representation.
    The binary vector is represented as a list of indices where the value is 1.
    """
    return [[1 if i in xi else 0 for i in range(n)] for xi in x]


def binary_distance(xi, n, xj=None):
    """
    Calculate the binary distance between all pairs of binary vectors in X.
    Expecting the binary vectors in "which" representation.
    The binary distance is defined as the sum of absolute differences
    between the corresponding elements of the vectors.
    """
    bin_xi = which2bin(xi, n=n)
    if xj is None:
        combs = list(itertools.combinations(bin_xi, 2))
        x_diff = [bin_d(pair[0], pair[1]) for pair in combs]
    else:
        bin_xj = which2bin(xj, n=n)
        x_diff = []
        for bxi in bin_xi:
            x_diff.append(min([bin_d(bxi, bxj) for bxj in bin_xj]))
    return x_diff


def set_matches(xi, xj):
    """
    Calculate the set matches between two sets of binary vectors.
    Expecting the binary vectors in "which" representation.
    The set distance is defined as the sum of exact matches
    """
    return len(set(which2dec(xi)) & set(which2dec(xj)))


def summarise(x):
    """
    Calculate the summary statistics for a list of numbers.
    The summary statistics include mean, standard deviation, median,
    minimum, and maximum of the binary distances.
    """
    return {
        "mean": np.mean(x),
        "std": np.std(x),
        "med": np.median(x),
        "min": np.min(x) if len(x) > 0 else np.NAN,
        "max": np.max(x) if len(x) > 0 else np.NAN,
    }


def mean_bin_dist(num_optima, n):
    full_list = list(itertools.product([0, 1], repeat=n))  # type: ignore
    xi = [np.where(x)[0] for x in full_list]
    xj = xi[:num_optima]
    binary_dist = binary_distance(xi, n=n, xj=xj)
    return summarise(binary_dist)


def get_run_data(n):
    api = wandb.Api()
    entity, project, group = "TAI_track", "hpo", "all"

    if group is not None:
        runs = api.runs(entity + "/" + project, filters={"group": group})
    else:
        runs = api.runs(entity + "/" + project)

    files = {}
    for run in runs:
        hpo_dict = defaultdict(dict)
        hpo_mins = defaultdict(dict)
        for key, value in run.summary.items():
            if "ensembles" in key:
                _, _, _, id, fold, _ = key.split("/")
                ensembles = json.loads(value)
                summary_dict = {}
                for k, v in ensembles.items():
                    if k == "1530":
                        # Dummy value
                        continue
                    X = [x for x, _ in v]
                    x_diff = binary_distance(X, n=n)
                    T = [t for _, t in v]
                    summary_dict[float(k)] = {
                        "number": len(v),
                        "bin_dist": summarise(x_diff),
                        "T": summarise(T),
                    }
                hpo_dict[id][fold] = summary_dict
                min_k = min(summary_dict.keys())
                hpo_mins[id][fold] = {
                    "X": [tuple(x) for x, _ in ensembles[str(min_k)]],
                    "T": [t for _, t in ensembles[str(min_k)]],
                    "min_k": min_k,
                }
        ensemble_file = f"ensembles_{run.id}.pkl"
        with open(ensemble_file, "wb") as f:
            pickle.dump(hpo_dict, f)
        min_file = f"ensembles_mins_{run.id}.pkl"
        with open(min_file, "wb") as f:
            pickle.dump(hpo_mins, f)
        files[run.id] = {
            "ensemble_file": ensemble_file,
            "min_file": min_file,
        }
    return files


parser = argparse.ArgumentParser()
parser.add_argument(
    "--run-id",
    type=str,
    default=None,
    help="The run id to fetch data for. If not provided, fetches data for all runs.",
)
parser.add_argument(
    "--n",
    type=int,
    default=9,
    help="The length of the binary vector. Default is 9.",
)


if __name__ == "__main__":
    args = parser.parse_args()
    if args.run_id is None:
        files = get_run_data(args.n)
    else:
        files = {
            args.run_id: {
                "ensemble_file": "ensembles_" + args.run_id + ".pkl",
                "min_file": "ensembles_mins_" + args.run_id + ".pkl",
            }
        }
    for run_id, file_paths in files.items():
        ensemble_file = file_paths["ensemble_file"]
        min_file = file_paths["min_file"]
        with open(ensemble_file, "rb") as f:
            ensembles = pickle.load(f)
        with open(min_file, "rb") as f:
            ensembles_mins = pickle.load(f)
        ensemble_mins = pickle.load(open(min_file, "rb"))
        ensemble_dict = pickle.load(open(ensemble_file, "rb"))
        results = []
        for id, fold_dicts in ensemble_mins.items():
            for fi in range(len(fold_dicts.keys())):
                fik = list(fold_dicts.keys())[fi]
                mins_fi = fold_dicts[fik]
                # Add stats for this fold
                fold_dict = {
                    "number": ensembles[id][fik][mins_fi["min_k"]]["number"],
                    "to": fik,
                    "id": id,
                    "from": fik,
                }
                fold_dict["bin_dist"] = mean_bin_dist(fold_dict["number"], n=args.n)
                results.append(pd.json_normalize(fold_dict))

                # Compute stats across folds
                for fj in range(fi + 1, len(fold_dicts.keys())):
                    fjk = list(fold_dicts.keys())[fj]
                    mins_fj = fold_dicts[fjk]
                    number = ensembles[id][fjk][mins_fj["min_k"]]["number"]
                    total_set_matches = set_matches(mins_fi["X"], mins_fj["X"])
                    results.append(
                        pd.json_normalize(
                            {
                                "number": number,
                                "from": fik,
                                "id": id,
                                "to": fjk,
                                "set_matches": total_set_matches,
                                "bin_dist": summarise(
                                    binary_distance(
                                        mins_fi["X"], n=args.n, xj=mins_fj["X"]
                                    )
                                ),
                            }
                        )
                    )
                    results.append(
                        pd.json_normalize(
                            {
                                "number": number,
                                "from": fjk,
                                "id": id,
                                "to": fik,
                                "set_matches": total_set_matches,
                                "bin_dist": summarise(
                                    binary_distance(
                                        mins_fj["X"], n=args.n, xj=mins_fi["X"]
                                    )
                                ),
                            }
                        )
                    )
        pd.concat(results).to_csv(f"{run_id}_results.csv")
