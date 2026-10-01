from typing import ClassVar

import numpy as np
from pymoo.core.individual import Individual

from mating_kernel.pymoo.parser.recording_parser import RecordingParser


class ReproductionParser(RecordingParser):
    record_args: ClassVar[list[str]] = ["selection", "survival", "crossover"]
    record_attributes: ClassVar[list[str]] = []
    record_retrieval: ClassVar[list[str]] = [
        "mating.selection",
        "mating.crossover",
        "survival",
    ]

    def __init__(self, ideal: np.ndarray | None = None):
        self.ideal = ideal

    def parse(self, data: dict[str, list]) -> list[dict]:
        # meta data about the recorded values
        meta = {}
        for k, v in data.items():
            meta[k] = v.pop(0)
            if k == "mating.selection":
                parents = [entry["output"] for entry in v]
            elif k == "mating.crossover":
                offspring = [entry["output"] for entry in v]
            elif (
                k == "survival"
            ):  # only need the last survival record, which is the final population
                survived = v[-1]["output"]
        assert len(parents) == len(offspring)
        ret_data = []
        for p, o in zip(parents, offspring):
            lineage = ReproductionParser.extract_lineage(
                p, o, n_offsprings=meta["mating.crossover"]["n_offsprings"]
            )
            for entry in lineage:
                # check if offspring survived in the final population
                entry_dict = ReproductionParser.parse_individual(
                    entry["offspring"], ideal=self.ideal
                )
                res_dict = {f"o_{k}": v for k, v in entry_dict.items()}
                for i, p in enumerate(entry["parents"]):
                    p_dict = ReproductionParser.parse_individual(p, ideal=self.ideal)
                    p_dict = {f"p{i}_{k}": v for k, v in p_dict.items()}
                    res_dict.update(p_dict)
                res_dict["survived"] = entry["offspring"] in survived
                ret_data.append(res_dict)

        return ret_data

    @staticmethod
    def parse_individual(ind: Individual, ideal: np.ndarray | None = None) -> dict:
        ind_dict = {"X": ind.X, "F": ind.F, **ind.data}
        if ideal is not None and len(ind.F) > 0:
            ind_dict["dist_to_ideal"] = np.linalg.norm(ind.F - ideal)
        return ind_dict

    @staticmethod
    def extract_lineage(parent_groups, offspring, n_offsprings):
        assert n_offsprings == 2, "Untested for n_offsprings != 2"
        n_matings = len(parent_groups)

        assert len(offspring) == n_matings * n_offsprings

        lineage = []

        for i, child in enumerate(offspring):
            offspring_slot = i // n_matings
            mating_idx = i % n_matings

            lineage.append(
                {
                    "offspring": child,
                    "offspring_slot": offspring_slot,
                    "mating": mating_idx,
                    "parents": list(parent_groups[mating_idx]),
                }
            )

        return lineage
