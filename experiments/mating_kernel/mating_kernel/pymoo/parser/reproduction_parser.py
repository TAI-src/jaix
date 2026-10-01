from mating_kernel.pymoo.parser.recording_parser import RecordingParser
from pymoo.core.individual import Individual


class ReproductionParser(RecordingParser):
    record_args = ["selection", "survival", "crossover"]
    record_attributes = []
    record_retrieval = ["mating.selection", "mating.crossover", "survival"]

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
                entry_dict = ReproductionParser.parse_individual(entry["offspring"])
                res_dict = {f"o_{k}": v for k, v in entry_dict.items()}
                for i, p in enumerate(entry["parents"]):
                    p_dict = ReproductionParser.parse_individual(p)
                    p_dict = {f"p{i}_{k}": v for k, v in p_dict.items()}
                    res_dict.update(p_dict)
                res_dict["survived"] = entry["offspring"] in survived
                ret_data.append(res_dict)

        return ret_data

    @staticmethod
    def parse_individual(ind: Individual) -> dict:
        return {"X": ind.X, "F": ind.F, **ind.data}

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
